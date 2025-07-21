from pathlib import Path

import numpy as np
import pandas as pd

from pycorn import PcUni6
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone

def get_series_from_data_dict(data_dictionary, target_key, data_key_list):
    try:
        # select the first injection as the injection timestamp
        inject_timestamp = data_dictionary[target_key]["Injection"]["data"][-1][0]
    except KeyError:
        inject_timestamp = 0

    data_series_list = []
    for data_key in data_key_list:
        data_array = np.array(data_dictionary[target_key][data_key]["data"]).astype(float)
        data_series = pd.Series(data=data_array[:, 1], index=data_array[:, 0])
        # remove duplicates
        data_series = data_series[~data_series.index.duplicated()]
        # offset by the infection_timestamp
        data_series.index -= inject_timestamp

        data_series_list.append(data_series)

    df = pd.concat(data_series_list, axis=1)
    df.columns = data_key_list
    return df


def import_xml_as_df(file_path: (str | Path), data_key_list: list = None, index: np.ndarray = None) -> pd.DataFrame:
    """
    Import the contents of a Unicorn Res/zip file into a pd.Dataframe

    Parameters
    ----------
    file_path : str or Path, Path to the res or zip file.
    data_key_list: list, optional, Keys to include in the DataFrame. Default: ["Cond", "UV", "Conc B"]
    index: np.ndarray, optional, Array of shape (1, ), to be used as index in the returned pd.DataFrame.
        Units are the same as the original data.

    Returns
    -------
    dataframe : pd.DataFrame

    """
    if data_key_list is None:
        data_key_list = ["Cond", "UV", "Conc B"]

    data_dictionary = PcUni6(file_path)
    data_dictionary.load_all_xml()

    target_key_list = [
        key for key in data_dictionary
        if "events" not in key
           and "Cond" in key
           and any(s.lower() in key.lower() for s in ["Tracer", "Injection", "Chrom", "Breakthrough", "Elution"])
           and any(["UV" in sub_key for sub_key in data_dictionary[key]])
    ]

    if any("breakthrough" in key.lower() for key in target_key_list):
        target_key_list = [key for key in target_key_list if "breakthrough" in key.lower()]

    if len(target_key_list) == 0:
        return None

    if "UV" not in data_dictionary[target_key_list[0]]:
        data_key_list.remove("UV")
        data_key_list.extend([sub_key for sub_key in data_dictionary[target_key_list[0]] if
                              ("UV" in sub_key and not "cell path" in sub_key)])

    series_list = [get_series_from_data_dict(data_dictionary, chrom, data_key_list) for chrom in target_key_list]

    if index is None:
        index = np.linspace(series_list[0].index.min(), series_list[0].index.max(), 100).round(3)

    # align and unify the index
    index = index[index < series_list[0].index.max()]

    if len(data_key_list) > 1:
        column_names = [(target_key, data_key) for target_key in target_key_list for data_key in data_key_list]
        column_names += [("temporary_insert", 0)]
    else:
        column_names = target_key_list
        column_names += ["temporary_insert"]

    # combine datapoints into one DataFrame and interpolate onto the unified index
    series_list.append(pd.Series(data=np.NaN, index=index))
    dataframe = pd.concat(series_list, axis=1, join="outer")
    dataframe.columns = column_names
    dataframe = dataframe.interpolate("index", )
    dataframe = dataframe.loc[index, column_names[:-1]]

    return dataframe

def get_metadata(data_dictionary: PcUni6|dict )-> dict[str, str]:
	"""
	Extracts metadata from the XML data of a Unicorn result file.
	must be called before get_chrom_logs because of the load_all_xm() method or any other function calling such method
	'creation_date' can be represented with .strftime('%d/%m/%y %H:%M:%S %Z') method

	Inputs:
		xml_data (PcUni6 or dict): Data structure containing loaded XML data from a Unicorn result file.

	Outputs:
		metadata (dict[str, str|datetime]): Dictionary containing extracted metadata fields:
			- 'result_name': Name of the result file
			- 'batch_ID': Batch ID associated with the result
			- 'creation_date': Creation date and time incl. timezone offset
			- 'system_name': Name of the system used for the run

	Notes:
		The function parses the XML content in xml_data['Result.xml'] and extracts relevant metadata fields.
		If a field is not found, its value will be None or 'not found'.
	"""
    
	# Parse the XML data from xml_data['Result.xml']
	xml_data = data_dictionary
	xml_data.load()
	root = ET.fromstring(xml_data['Result.xml'])


	# Find name
	name_elem = root.find('.//Name')
	result_name = name_elem.text if name_elem is not None else None

	# Find batch-id
	batch_id_elem = root.find('.//BatchId')
	batch_id = batch_id_elem.text if batch_id_elem is not None else None


	# Find creation date
	created_elem = root.find('.//Created')
	created_timestamp = created_elem.text if created_elem is not None else None
	
	utc_offset_elem = root.find('.//CreatedUtcOffsetMinutes')
	utc_offset_minutes = int(utc_offset_elem.text) if utc_offset_elem is not None else 0

	if created_timestamp:
		dt = datetime.fromisoformat(created_timestamp)
		# Apply UTC offset
		dt = dt.replace(tzinfo=timezone(timedelta(minutes=utc_offset_minutes)))
	else:
		dt = 'not found'



	# Find System Name
	system_name_elem = root.find('.//SystemName')
	system_name = system_name_elem.text if system_name_elem is not None else None

	col_name = xml_data['ColumnTypeData']['Xml']['ColumnTypes']['ColumnType']['Name']

	# print("Name: \t\t", result_name)
	# pritn("Column: \t", col_name) 
	# print("Created:\t", dt.strftime('%d/%m/%y %H:%M:%S %Z'))
	# print("Batch ID:\t", batch_id)
	# print("SystemName:\t", system_name)


	metadata = {
		'result_name': result_name,
		'column': col_name,
		'creation_date': dt,
		'batch_ID': batch_id,
		'system_name': system_name
	}
	return metadata

def get_chrom_from_data_dict(data_dictionary: PcUni6|dict, chromatogram_str, traces_list) -> pd.DataFrame:
    """""
    extract the chormatogram data from a data_dictionary for a given chromatogram_str and traces_list
    works one chromatogram at a time
    
    Inputs:
		data_dictionary (dict): dictionary containing Unicorn results with all chromatograms, as prepared by PcUni6()
		chromatogram_str (str): string representing the chromatogram name to extract data from
		traces_list (list): list of traces to extract from the chromatogram
        
    Outputs:
		df (pd.DataFrame): DataFrame containing the data for the given chromatogram and traces
        
	"""""
	
    inject_timestamp = 0
    
    traces_not_in_chromatogram: list[str] = [] # needed to handle case where some traces are only present in some chormatograms
    data_series_list = []

    for data_key in traces_list:
        if not(any(key.lower() == data_key.lower() for key in data_dictionary[chromatogram_str].keys())):
            traces_not_in_chromatogram.append(data_key)
            continue
        
        data_array = np.array(data_dictionary[chromatogram_str][data_key]["data"])
        if data_array.size == 0:
            x_data, y_data = np.array([np.nan]), np.array([np.nan])
        else: 
            x_data = data_array[:, 0].astype(float)
            y_data = data_array[:, 1]
        # data_series = pd.Series(data=y_data, index=x_data)
        try:
            data_series = pd.Series(data=y_data, index=x_data, dtype=float)
        except ValueError:
            data_series = pd.Series(data=y_data, index=x_data, dtype=pd.StringDtype())
        # remove duplicates
        data_series = data_series[~data_series.index.duplicated()]
        # offset by the injection_timestamp
        data_series.index -= inject_timestamp

        # add 'True' to the Injection data as mark should there be multiple injections
        if data_key == "Injection":
            data_series.replace(np.nan, True, inplace=True)

        data_series_list.append(data_series)

    try: # case where only one trace (besides "Injection") is present
        df = pd.concat(data_series_list, axis=1)
    except ValueError:
        pass
    
	# removes all all elems of 'traces_not_in_chromatogram' from 'traces_list' in case-insensitive manner
    traces_list_new = [t for t in traces_list if not any(t.lower() == rem.lower() for rem in traces_not_in_chromatogram)]
    df.columns = traces_list_new
    return df

def get_chrom_logs(data_dictionary: PcUni6|dict, **kwargs) -> tuple[pd.DataFrame, pd.DataFrame]:
	"""
	Extracts chromatogram and log data from a dictionary containing Unicorn results with all chromatograms, as prepared by PcUni6(), and returns two DataFrames: one for chromatograms and one for logs
	if no 'Injection' is provided in 'traces' (**kwargs), chromatogram will not be adjusted to the first injection, if no 'traces' is provided it will depended on if any injection was done

	Inputs:
		data_dictionary (dict): dictionary containing Unicorn results with all chromatograms, as prepared by PcUni6()
		kwargs (dict)
		chromatograms (list[int]): list of chromatogram names to import. If not provided, all chromatograms will be used.
		traces (list[str]): list of traces to import. If not provided, all traces will be used.
		interpolate_threshold (int): to control interpolation behavior, chromatograms with less rows than this will not be interpolated.
	
	Output:
		chromatogram_df (pd.DataFrame): DataFrame containing all chromatograms with aligned data
		log_df (pd.DataFrame):  DataFrame containing log data
	"""

	# Load all to data_dictionary
	data_dictionary.load_all_xml()

	# if 'chromatograms' list is not provided, all chromatograms from the data_dictionary will be used
	chromatograms: list[int] = kwargs.get('chromatograms', [key for key in list(data_dictionary.keys()) if not key.lower().endswith('.xml_dict')])
	# Threshold for interpolation, can be adjusted, chromatograms with less than this number of rows will not be interpolated
	interpolate_threshold: int = kwargs.get('interpolate_threshold',10)

	# Initialize empty DataFrames and lists to store chromatogram and log data
	chromatogram_df: pd.DataFrame = pd.DataFrame()
	log_df: pd.DataFrame = pd.DataFrame()
	chrom_dfs_list: list = []
	log_dfs_list: list = []
		
	for chrom_name in chromatograms:
		print(f"processing \t{chrom_name}")
		traces = kwargs.get('traces', list(data_dictionary[chrom_name].keys()))
		 
		df = get_chrom_from_data_dict(data_dictionary, chrom_name, traces).sort_index().dropna(how='all')

		df_logs = df.select_dtypes(exclude=['float64'])
		df_logs.dropna(how='all', inplace=True)
		

		df = df.select_dtypes(include=['float64'])

		# interpolate to allign all y-values (chromatograms) to the same x-values (mL coordinates)
		if df.shape[0] > interpolate_threshold:
			df.count(axis=0).median()  # Count non-NA/null observations in each column

			# Generate new index with median df index length, preserving min and max of original df index
			new_index = np.linspace(df.index.min(), df.index.max(), int(df.count(axis=0).median()))

			# Interpolate all columns to the new index (merge two indexes, interpolate the resuling missing values, select only the new index)
			df_interpolated = df.reindex(df.index.union(new_index)).interpolate(method='index').loc[new_index]

			# Set the index name to match the original if needed
			df_interpolated.index.name = df.index.name
		else:
			df_interpolated = df.copy()

		# append the DataFrames to the respective lists for later concatenation
		log_dfs_list.append(df_logs)
		chrom_dfs_list.append(df_interpolated)
		

	#merge all chromatograms and logs into respective DataFrames
	chromatogram_df = pd.concat(chrom_dfs_list, axis=0)
	log_df = pd.concat(log_dfs_list, axis=0)

	if "Injection" in traces:
	# find first injection and reajust the index
		first_injection_idx = log_df['Injection'][log_df['Injection'] == True].index.min()
		chromatogram_df.index = chromatogram_df.index - first_injection_idx
		log_df.index = log_df.index - first_injection_idx
	
	chromatogram_df.sort_index(inplace=True)
	log_df.sort_index(inplace=True)

	return chromatogram_df, log_df


