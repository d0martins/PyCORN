#standard library
from collections.abc import Sequence
from pathlib import Path
from datetime import datetime, timedelta, timezone
from math import pi

# third-party
from pycorn.utils import get_between_logs
from pycorn import PcUni6
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

# local
from aktachromatogram.model_dataclass import Result, ResultBatch



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
	must be called before get_chrom (because of the load_all_xm() method) or before any other function calling such method
	'creation_date' can be represented with .strftime('%d/%m/%y %H:%M:%S %Z') method

	Inputs:
		xml_data (PcUni6 or dict): Data structure containing loaded XML data from a Unicorn result file.

	Outputs:
		metadata (dict[str, str|datetime]): Dictionary containing extracted metadata fields:
			- 'result_name': Name of the result file
			- 'column': column name as saved in Unicorn column list
			- 'column_bed_height_cm': as saved in Unicorn column list
 			- 'column_diameter_cm': as saved in Unicorn column list
 			- 'column_volume_mL': as saved in Unicorn column list
			- 'creation_date': Creation date and time incl. timezone offset
			- 'batch_ID': Batch ID associated with the result
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
	
	# Find column data if available
	if (xml_data['ColumnTypeData']['Xml'] is not None) and ("ColumnType" in xml_data['ColumnTypeData']["Xml"]["ColumnTypes"].keys()):
		# for runs started in manual mode the xml_data['ColumnTypeData']['Xml'] is empyt (None)
		# for runs started by a method the xml_data['ColumnTypeData']['Xml'] exists with the xml-schema only and
		# nothing else, therefore the additional test is needed for method runs without column information
		col_name = xml_data['ColumnTypeData']['Xml']['ColumnTypes']['ColumnType']['Name']
		col_bed_height_unit = str(xml_data['ColumnTypeData']["Xml"]["ColumnTypes"]["ColumnType"]["BedHeightUnit"])
		col_bed_height_key = "column_bed_height_" + col_bed_height_unit
		col_bed_height = float(xml_data['ColumnTypeData']["Xml"]["ColumnTypes"]["ColumnType"]["BedHeight"])
		col_diameter_unit = str(xml_data['ColumnTypeData']["Xml"]["ColumnTypes"]["ColumnType"]["Hardware"]["DiameterUnit"])
		col_diameter_key = "column_diamter_" + col_diameter_unit
		col_diameter = float(xml_data['ColumnTypeData']["Xml"]["ColumnTypes"]["ColumnType"]["Hardware"]["Diameter"])
		if col_bed_height_unit.lower() == "cm" and col_diameter_unit.lower() == "cm":
			col_volume_key = "column_volume_mL"
		else:
			col_volume_key = "column_volume_" + col_bed_height_unit + "*" +  col_diameter_unit + "^2"
		col_volume = round(pi * (col_diameter/2)**2 * col_bed_height, 4)
	else:
		col_name = None
		col_bed_height_key = "column_bed_height"
		col_bed_height = None
		col_diameter_key = "column_diameter"
		col_diameter = None
		col_volume_key = "column_volume"
		col_volume = None

	# print("Name: \t\t", result_name)
	# pritn("Column: \t", col_name) 
	# print("Created:\t", dt.strftime('%d/%m/%y %H:%M:%S %Z'))
	# print("Batch ID:\t", batch_id)
	# print("SystemName:\t", system_name)


	metadata = {
		'result_name': result_name,
		'column': col_name,
		col_bed_height_key: col_bed_height,
		col_diameter_key: col_diameter,
		col_volume_key: col_volume,
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

    try: 
        df = pd.concat(data_series_list, axis=1)
    except ValueError: # case where only one trace (besides "Injection") is present
        pass
    
	# removes all all elems of 'traces_not_in_chromatogram' from 'traces_list' in case-insensitive manner
    traces_list_new = [t for t in traces_list if not any(t.lower() == rem.lower() for rem in traces_not_in_chromatogram)]
    df.columns = traces_list_new
    return df


def get_chrom(data_dictionary: PcUni6|dict, reduce_interpolate: bool = False, **kwargs) -> tuple[pd.DataFrame, pd.DataFrame]:
	"""
	Extracts chromatogram from a dictionary containing Unicorn results with multiple chromatograms, as prepared by PcUni6(), and returns two DataFrame: one with consolidated data from chromatogram(s) and one with the logs of fractions and injection(s)
	if no 'Injection' is provided in 'traces' (**kwargs), chromatogram will not be adjusted to the first injection, if no 'traces' is provided it will depended on if a 'Injection' trace exist (i.e. injection was done)

	Inputs:
		data_dictionary (dict): dictionary containing Unicorn results with all chromatograms, as prepared by PcUni6()
		reduce_interpolate (bool): False by default, whether to reduce/interpolate the data to a frequency equal to the mean of the trace frequency (see note below)
		kwargs (dict)
		chromatograms (list[int]): list of chromatogram names to import. If not provided, all chromatograms will be used.
		traces (list[str]): list of traces to import. If not provided, all traces will be used.
		interpolate_threshold (int): to control interpolation behavior, chromatograms with less rows than this will not be interpolated.
        which_injection (int): which injection mark to use (0-indexed injection mark, -1 for last), default is first injeciton mark
	
	Output:
		chromatogram_df (pd.DataFrame): DataFrame containing all chromatograms with aligned data
		frac_log_df (pd.DataFrame):  DataFrame containing the fractions (and injection(s)) logs
	"""

	# Load all to data_dictionary
	data_dictionary.load_all_xml()

	# if 'chromatograms' list is not provided, all chromatograms from the data_dictionary will be used
	chromatograms: list[int] = kwargs.get('chromatograms', [key for key in list(data_dictionary.keys()) if not key.lower().endswith('.xml_dict')])
	# Threshold for interpolation, can be adjusted, chromatograms with less than this number of rows will not be interpolated
	interpolate_threshold: int = kwargs.get('interpolate_threshold', 10)

	# Initialize empty DataFrames and lists to store chromatogram and log data
	chromatogram_df: pd.DataFrame = pd.DataFrame()
	log_df: pd.DataFrame = pd.DataFrame()
	chrom_dfs_list: list = []
	log_dfs_list: list = []
		
	for chrom_name in chromatograms:
		#print(f"processing \t{chrom_name}")
		traces = kwargs.get('traces', list(data_dictionary[chrom_name].keys()))
		traces = [trace for trace in traces if trace.lower() not in "Run Log".lower()]
		 
		df = get_chrom_from_data_dict(data_dictionary, chrom_name, traces).sort_index().dropna(how='all')

		df_logs = df.select_dtypes(exclude=['float64'])
		df_logs.dropna(how='all', inplace=True)
		

		df = df.select_dtypes(include=['float64'])

		if reduce_interpolate: # interpolate to allign all y-values (chromatograms) to the same x-values (mL coordinates)
			if df.shape[0] > interpolate_threshold:
				df.count(axis=0).median()  # Count non-NA/null observations in each column

				# Generate new index with median df index length, preserving min and max of original df index
				new_index = np.linspace(df.index.min(), df.index.max(), int(df.count(axis=0).median()))
				
				# Interpolate all columns to the new index (merge two indexes, interpolate the resuling missing values, select only the new index)
				df_interpolated = df.reindex(df.index.union(new_index)).interpolate(method='index').loc[new_index]

				# Set the index name to match the original if needed
				df_interpolated.index.name = df.index.name
		
		df_interpolated = df.copy()

		# append the DataFrames to the respective lists for later concatenation
		log_dfs_list.append(df_logs)
		chrom_dfs_list.append(df_interpolated)
	chromatogram_df = pd.concat(chrom_dfs_list, axis=0) #merge all chromatograms and logs into respective DataFrames
	log_df = pd.concat(log_dfs_list, axis=0)
	
	if "Injection" in traces:
		injection_selected = int(kwargs.get('which_injection', 0)) # find first injection and reajust the index    
		injection_idx = log_df['Injection'][log_df['Injection'] == True].index[injection_selected]
		chromatogram_df.index = chromatogram_df.index - injection_idx
		log_df.index = log_df.index - injection_idx
	
	chromatogram_df.sort_index(inplace=True)
	frac_log_df = log_df.sort_index().copy()
	return chromatogram_df, frac_log_df


def get_full_log(data_dictionary: PcUni6|dict, **kwargs) -> pd.DataFrame:
	"""
	extract the full log
	
	Inputs:
		data_dictionary (dict): dictionary containing Unicorn results with all chromatograms, as prepared by PcUni6()
        kwargs (dict)
        which_injection (int): which injection mark to use (0-indexed injection mark, -1 for last), default is first injeciton mark
	Outputs:
		full_log (pd.DataFrame): df with the full log

	"""
	chromatograms: list[int] = [key for key in list(data_dictionary.keys()) if key.lower().endswith('.xml_dict')]
	frames = []
	for chrom in chromatograms:
		xml_chrom = data_dictionary[chrom]
		events = xml_chrom["Chromatogram"]["EventCurves"]["EventCurve"]
		if isinstance(events, dict):
			
			if events["@EventCurveType"] == "Logbook":
				log = events["Events"]["Event"]
				partial_log = pd.DataFrame.from_records(log)
		elif isinstance(events, list):
			for event_curve in events:
				if event_curve['@EventCurveType'] == "Logbook":
					log = event_curve["Events"]["Event"]
					partial_log = pd.DataFrame.from_records(log)
		else:
			raise TypeError(f"'dict'  or 'list' type expected at xml_data['{chrom}']['Chromatogram']['EventCurves']['EventCurve']")
		
		frames.append(partial_log\
				.astype({"EventTime": "float", "EventVolume": "float"}))
		
		# get injection points, try to find in the current chromatogram if exists
		try:
			for event in data_dictionary[chrom]["Chromatogram"]["EventCurves"]["EventCurve"]:
				if event["@EventCurveType"] == "Injection":
					event_entry = event["Events"]["Event"]
					
					if isinstance(event_entry, dict):
						injection = {
							"volume_ml" : float(event_entry["EventVolume"]),
							"time_min" : float(event_entry["EventTime"])}
					elif isinstance(event_entry, list):
						# selects the first injection done
						injection_selected = int(kwargs.get("which_injection", 0))
						injection = {
							"volume_ml" : float(event_entry[injection_selected]["EventVolume"]),
							"time_min" : float(event_entry[injection_selected]["EventTime"])}
		except TypeError:
			# print(f"no injectin found in {chrom}")
			pass

	full_log = pd.concat(frames)
	full_log.columns = [col.replace("@", "") for col in full_log]
	full_log.sort_values(by="EventTime", inplace=True, ignore_index=True)
	full_log['EventFullText'] = full_log[['EventType', 'EventSubType', 'EventText']].apply(lambda row: ' '.join(row.astype(str)), axis=1)
	full_log = full_log[["EventVolume", "EventTime", "EventType", "EventSubType", "EventText", "InstructionFeedback", "EventFullText"]]

	# IF no injection mark was found, 0 (zero) will be set by default,
	# in order be able to adjust "EventTime" and "EventVolume" in any case (with and without injectio mark)
	try:
		assert(injection)
	except NameError:
		injection = {
			"volume_ml" : float(0),
			"time_min" : float(0)}

	full_log["EventTime"] = full_log["EventTime"] - injection["time_min"]
	full_log["EventVolume"] = full_log["EventVolume"] - injection["volume_ml"]

	return full_log


def get_frac_vol(log_df: pd.DataFrame)-> pd.DataFrame:
	"""
	calculate fraction volume from 'log_df' (short version)
	
	Inputs:
		log_df (pd.DataFrame): df with fraction names and volume coordinates as prepared by _, log_df=get_chrom_logs() from pycorn.utils
	
	Outputs:
		frac_vol_df (pd.DataFrame): df with only fractions, their volume and start point (=log_df.index)
	"""

	df = log_df.copy()

	df.dropna(subset = "Fractions", inplace= True)
	df["frac_start_volume_ml"] = df.index
	df["fraction_volume"] = - df["frac_start_volume_ml"].diff(periods = -1) # "Fractions" marks the beginning of the fraction therefore negative and period = -1 (difference to following line) is needed

	df.drop(df[df["Fractions"].str.contains("Waste|Frac", case=False)].index, inplace=True) # drops lines with 'frac' or 'waste' as fraction name

	frac_vol_df = df[["frac_start_volume_ml", "Fractions", "fraction_volume"]]

	return frac_vol_df


def get_between_logs(full_log_df: pd.DataFrame, start_end_text: list[str], lookup_col: str = "EventFullText", return_col: str = "EventVolume"):
    """
    Extract values from a DataFrame corresponding to two log entries matching given text patterns.

    recommended usage: <some_df>.loc[slice(*edges)] (volume coordinate expected in some_df.index)
    
    Parameters
    ----------
    full_log_df : pd.DataFrame
        Log data containing at least the columns specified by `lookup_col` and `return_col`.
    start_end_text : list of str
        List of two string patterns (start and end markers) to locate within `lookup_col`.
        The first match corresponds to the start event, and the second to the end event.
    lookup_col : str, optional
        Name of the column in which to search for `start_end_text` patterns.
        Default is "EventFullText", alternatives are `EventText`, `EventType`, 
        `EventSubType` and `InstructionFeedback`
    return_col : str, optional
        Name of the column from which to extract values corresponding to matched rows.
        Default is "EventVolume"; alternative is `EventTime`

    Returns
    -------
    edges : list of float
        List containing the two extracted values (start and end) from `return_col`.

    Raises
    ------
    IndexError
        If no matching rows are found for either pattern.
    KeyError
        If `lookup_col` or `return_col` do not exist in `full_log_df`.

    Notes
    -----
    This function searches case-insensitively for each string pattern in `start_end_text`
    within the specified column, retrieves the first match of each, and returns their
    associated numeric values. Useful for isolating chromatographic or process segments
    between two event markers in log data.
    """
    edges = list(map(lambda pattern_str: full_log_df[full_log_df[lookup_col].str.contains(pattern_str, case=False)][return_col].values[0], start_end_text))
    return edges


def get_fracs_between_logs(full_log_df: pd.DataFrame, frac_df: pd.DataFrame, start_end_text:list[str]):
	"""
	gets first and last fraction between two log events ("EventFullText" column)
	
	Parameters
	---------
	full_log_df : pd.DataFrame
		run log; must contain columns ["EventFullText", "EventVolume"]
	frac_df : pd.DataFrame
		fraction table; must contain columns ["frac_start_volume_ml"]
	
	Returns
	-------
	log_frac_series : pd.Series
		series with all fractions between the text in the logs (index is injection-zeroed elution volume (mL))

	"""
	edges = get_between_logs(full_log_df, start_end_text)
	print(edges)
	closest_indices = [(frac_df["frac_start_volume_ml"] - target).abs().idxmin() for target in edges] # computes the absolute distance between each row and the target, gets the index of the smallest distance (i.e., the closest match)
	print(closest_indices)
	log_frac_series = frac_df.loc[closest_indices[0]:closest_indices[1], "Fractions"].str.replace(".", "")
	return log_frac_series


def get_start_end_frac(dfm: pd.DataFrame, chrom_data_dict: dict, expid: str, params_exp: dict, inbetween_fractions = True, **kwargs):
	"""
	get the fractions between specific log event text which are ALSO present in the lims data ('dfm')

	Parameters
	----------
	dfm : pd.DataFrame
		df (multiindex) with LIMS data
	expid : str
		experiment ID/key, e.g. "250709KLSE_A"
	params_exp : dict
		keys include "cycle", "resin" "FcXP", "vol_col_mL", "flow_rate_mLmin-1", "fraction_size_bt_mL"
	inbetween_fractions : bool, default True
		whether to report fractions between fist and last fraction or not
	**kwargs
		start_end_text (list[str]): 2-element list with strings for matching in full_log["EventFullText"] column marking start and end of desired portion
		dest_key (list[str]): 2-element list with destination keys to be updated on the 'params_exp' dict

	Returns
	-------
		start_end_frac : pd.Series
		series with start- and end-fraction name (optionally also fractions inbetween)
	"""
	start_end_text = kwargs.get("start_end_text", ["Method BlockStart Block Sample Inlet",  "Method BlockStart Phase Wash after Load"])

	# cycle = params_exp["cycle"]

	# LIMS fractions, query of
	_, df_above_loq = loq_splitter(dfm)
	lims_frac_series = df_above_loq[df_above_loq["expID"]==expid]["sample"]

	# Fractions in the log (from Unicorn result file), query of
	full_log_df = chrom_data_dict["full_log"]
	frac_df = chrom_data_dict["frac_vol"]
	log_frac_series = get_fracs_between_logs(full_log_df, frac_df, start_end_text)

	# intersction of Log Phase log (Unicorn file) ∩ LIMS
	union_frac_series = log_frac_series[log_frac_series.isin(lims_frac_series)]

	if inbetween_fractions:
		start_end_frac = union_frac_series
	else: 
		start_end_frac = union_frac_series.iloc[[0,-1]]

	return start_end_frac


def extract_trace_series(
    *,
    batch: ResultBatch | None = None,
    results: Sequence[Result] | None = None,
    trace_key: str = "DeltaC pressure",
    start_end_text: Sequence[str] | None = None,
    top_cycle_cutoff: int | None = None
) -> pd.Series:
    """
    Extract a trace (e.g., pressure, UV, ...) from multiple chromatography cycles
    and return as a Series indexed by cycle_count.
    
    Each value in the Series contains the full trace data (as pd.Series) for that cycle.

    Parameters
    ----------
    batch : ResultBatch, optional
        A ResultBatch instance with sorted results. Mutually exclusive with `results`.

    results : list[Result], optional
        A filtered list of Result instances. Mutually exclusive with `batch`.

    trace_key : str, default="DeltaC pressure"
        The chromatogram trace to extract.

    start_end_text : list[str], optional
        Custom start and end event markers to extract trace segment.
        If None, uses the full chromatogram.

    top_cycle_cutoff : int, optional
        Max number of cycles to process. Defaults to all.

    Returns
    -------
    pd.Series
        Series indexed by cycle_count, where each value is a pd.Series containing
        the trace data for that cycle.
        
    Examples
    --------
    >>> traces = extract_trace_series(
    ...     batch=batch,
    ...     trace_key="DeltaC pressure",
    ...     start_end_text=["BlockStart Block Direct sample injection_1", "(Column Wash)"]
    ... )
    >>> # Access trace for cycle 5
    >>> cycle_5_pressure = traces[5]
    >>> # Get median for cycle 5
    >>> median_5 = traces[5].median()
    """
    if (batch is None and results is None) or (batch and results):
        raise ValueError("Provide either 'batch' or 'results', but not both.")

    if results is None:
        results = batch.results

    if top_cycle_cutoff:
        results = results[:top_cycle_cutoff]

    trace_data = {}

    for result in results:
        try:
            chrom_full = result.chrom
            full_log = result.full_log

            if start_end_text is not None:
                edges = get_between_logs(full_log, start_end_text)
                chrom = chrom_full.loc[slice(*edges)]
            else:
                chrom = chrom_full

            df_cycle = chrom.loc[:, trace_key].dropna()
            
            if len(df_cycle) > 0:
                trace_data[result.cycle_count] = df_cycle

        except Exception as e:
            print(f"[{result.cycle_count}] Error extracting trace: {e}")
            continue

    return pd.Series(trace_data, name=trace_key)



def interpolate_to_column(
    df: pd.DataFrame,
    source_col: str,
    target_col: str,
    new_col: str | None = None,
    method: str = 'linear',
    fill_value: float | str = np.nan
) -> pd.Series | pd.DataFrame:
    """
    Interpolate source_col values to align with the non-NaN indices of target_col.
    
    This is useful when two measurement columns have data at different time points
    (misaligned indices) and you need to align them for mathematical operations.
    
    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe with misaligned columns.
    source_col : str
        The column to interpolate (e.g., "Sample flow").
    target_col : str
        The column whose non-NaN indices define the output alignment (e.g., "DeltaC pressure").
    new_col : str, optional
        If provided, adds the interpolated values as a new column to the dataframe
        and returns the modified dataframe. If None, returns only the interpolated Series.
    method : str, default='linear'
        Interpolation method. Options: 'linear', 'nearest', 'cubic', etc.
    fill_value : float or str, default=np.nan
        How to handle extrapolation beyond the source data range.
        - np.nan: use NaN for out-of-bounds values
        - 'extrapolate': extend the interpolation beyond bounds
        - float: use a specific fill value
    
    Returns
    -------
    pd.Series or pd.DataFrame
        If new_col is None: returns interpolated Series aligned to target_col's non-NaN index.
        If new_col is provided: returns DataFrame with new column added.
    
    Examples
    --------
    >>> # Get interpolated series
    >>> flow_aligned = interpolate_to_column(df, "Sample flow", "DeltaC pressure")
    >>> 
    >>> # Add as new column
    >>> df = interpolate_to_column(df, "Sample flow", "DeltaC pressure", new_col="Sample flow (aligned)")
    >>> 
    >>> # Use for calculations
    >>> flow_interp = interpolate_to_column(df, "Sample flow", "DeltaC pressure")
    >>> df["pressure_per_flow"] = df["DeltaC pressure"].dropna() / flow_interp
    """
    # Get non-NaN indices from target column
    target_index = df[target_col].dropna().index.values
    
    # Get source data (drop NaN)
    source_series = df[source_col].dropna()
    source_index = source_series.index.values
    source_values = source_series.values
    
    # Interpolate source values to target indices
    if method == 'linear':
        # Handle fill_value for np.interp
        if fill_value == 'extrapolate':
            # np.interp extrapolates by default at boundaries
            interpolated = np.interp(target_index, source_index, source_values)
        else:
            interpolated = np.interp(target_index, source_index, source_values)
            # Mask values outside source range
            out_of_bounds = (target_index < source_index.min()) | (target_index > source_index.max())
            if isinstance(fill_value, (int, float)):
                interpolated[out_of_bounds] = fill_value
            else:  # np.nan or other
                interpolated[out_of_bounds] = np.nan
    else:
        # Use scipy for other interpolation methods
        from scipy.interpolate import interp1d
        f = interp1d(source_index, source_values, kind=method, 
                     bounds_error=False, fill_value=fill_value)
        interpolated = f(target_index)
    
    # Create output series
    result_series = pd.Series(interpolated, index=target_index, name=source_col)
    
    if new_col is not None:
        # Add as new column to dataframe
        df = df.copy()
        df[new_col] = np.nan
        df.loc[target_index, new_col] = interpolated
        return df
    else:
        return result_series
