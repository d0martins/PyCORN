import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

import pandas as pd
from scipy import fft, interpolate, integrate

from scipy.stats import exponnorm

import warnings


def baseline_corrector(x: list[float]|np.ndarray, y: list[float]|np.ndarray, baseline_method: str = "spline", **kwargs: dict) -> list[float]|np.ndarray:
	"""
	calculate the baseline for a given peak (also clip small negative value to zero).
	Anchor points are selected from the derivative of y-data below a give threshold.
	The baseline is calculated by one of three methods: spline (default), interpolation or percentil
	Spline works well for a long and stable baselines. 
	Interpolation is backup option for short chromatograms and/or if the signal does not reach the initial baseline pre-peak
	Percentil takes a given level as the baseline (default is 5% percentil)

	Inputs:
		x (list[float]|np.ndarray): x-coordinate data e.g. eluted volume
		y (list[float]|np.ndarray): y-coordinate data e.g. Conductivity signal
		method (str): baseline calculation method: spline, interpolation or percentil
		**kwargs (dict):
		percentil_level (int): percentil level to take as baseline

	Output:
		y_corrected (list[float]|np.ndarray): baseline-corrected y-coordinate data e.g. baseline corrected signal

	NOTES:
		THRESHOLD to select anchor points
		np.average(dy) preferred over np.median(dy) because
			(1) 'dy' ofter initially overshoots then undershoots in the peak front and tail, respecitively, thus median (middle of population) for non-gaussian (symmetric) peaks
			(2) the 'dy' takes both positive and negative values and therefore the median threshold could lead to the inclusion of some part of the peak as anchor points
			(3) median(abs(dy)) is also not possible, since an anchor point is ofter identified at peak apex
	"""
	# derivative of y-data
	window_length: int = kwargs.get("anchors_window_length", 11)
	polyorder: int = kwargs.get("anchors_polyorder", 3)
	dy = np.gradient(savgol_filter(y, window_length, polyorder))

	#THRESHOLD (see note in docstring above)
	threshold = abs(np.average(dy)) 			#1e-3 as alternative value

	anchors = np.abs(dy) < threshold 			# |dy| small → baseline

	if baseline_method.lower() == "spline":
		# SPLINE BASELINE, derivative-based anchor-points
		# default option, works with long a stable baseline
		spline = UnivariateSpline(x[anchors], y[anchors], s=0) 	# smooth spline
		baseline = spline(x)
	elif baseline_method.lower() == "interpolation":
		# STRAIGHT LINE BASELINE, derivative-based anchor-points
		# backup for case without long and stable baseline (either short chromatogram or signal does not reach baseline)
		baseline = np.interp(x = x, xp = x[anchors], fp = y[anchors])
	elif baseline_method.lower() == "percentil":
		# simple constant subtraction and reset it as baseline
		percentil = kwargs.get("percentil_level", 5)
		baseline = np.percentile(y, percentil)
	else:
		raise ValueError("unexpected base line calculation method; 'interpolation', 'spline' or 'percentil' is expected!")

	y_corrected = y - baseline
	y_corrected[y_corrected < 0] = 0

	## plotting for visual check
	## 
	# plt.plot(x, baseline, '--', label='baseline')
	# plt.plot(x, y, label = "true data")
	# plt.scatter(x[anchors], y[anchors], label = "anchors")
	# plt.xlabel("x-data")
	# plt.ylabel("y-data")
	# plt.legend()
	# plt.show()
	# plt.close()
	# plt.plot(x, dy, label = "dy")
	# plt.axhline(0, color = "black", label = "zero-line")
	# plt.xlabel("x-data")
	# plt.ylabel("derivative of y-data")
	# plt.legend()
	# plt.show()
	# plt.close()

	return y_corrected, baseline


def sampler_dist_builder(x: (list[float]|np.ndarray), y: (list[float]|np.ndarray), plotting: bool = False) -> list|np.ndarray:
    """
    
    """
	# reconstitution to 1-D
    dx  = np.gradient(x)										# bin widths
    pdf = y.copy()
    pdf /= np.sum(pdf * dx)								# normalise (defensive)
    #print(np.sum(pdf * dx)) 							# check, integral of pdf must be 1

	# build CDF
    cdf = np.cumsum(pdf * dx)
    #print(cdf.max())											# check, max/last must be 1
    
	# resampling
    n:int = 100_000
    u: np.ndarray = np.random.rand(n)
    sample = np.interp(u, cdf, x)                 		# 1-D array
    

    if plotting: ## plotting for visual check
        plt.hist(sample, bins=60, density=True, alpha=0.4, label='resampled binned data')
        scaling_factor: float = pdf.max()/y_corrected.max()
        plt.plot(x, y*scaling_factor, label = "scaled true data")
        plt.plot(x, pdf, lw=2, label='reconstituted PDF', linestyle = "dotted")
        plt.legend(); plt.xlabel('x'); plt.ylabel('density'); plt.show()

    return sample, pdf


def emg_fitter(x: (list[float]|np.ndarray), y: (list[float]|np.ndarray), **kwargs)-> tuple[dict[str, float], dict[int, float], dict[str, float]]:
	"""
	fit exponentially-modified gaussian from xy-data and returns dictionaries with each central moments, raw moments and distribution parameters
	Plotting option (defaut Off) availabe for visualization/troubleshooting.

	Inputs:
		x (list[float]|np.ndarray): x-coordinate data e.g. eluted volume
		y (list[float]|np.ndarray): y-coordinate data e.g. Conductivity signal
		kwargs:
		plotting (bool): whether to plot (or not) the peak and EMG dist (PDF)
		baseline_method (str): baseline calculation method can be either 'spline' or 'interpolation'. 
		
	Output:
		dist_moments (dict[str, float]): dict with central moments (up to fourth) 
		dist_moments_raw (dict[int, float]): dict with raw (non-central) moments
		dist_params (dict[str, float]): dict with distribution parameters


	NOTES:
		NOMENCLATURE
		scipy.stats.exponnorm used a sligtly different parameter nomenclature than other websites/books
		(SciPy docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.exponnorm.html)
		loc, scale and K from the SciPy docs translates to mu (μ), sigma (σ), and k with k = 1/(sigma*lam) or lam = 1/(k*sigma) in this script and in the Wikipedia EMG page (lam being the exponential rate λ)
	"""

	plotting: bool = kwargs.get("plotting", False ) # whether to plot or not
	baseline_method: str = kwargs.get("baseline_method", "spline") # 'spline' as default methhod for baseline finding/correction

	y_corrected, baseline = baseline_corrector(x, y, baseline_method = baseline_method)
	sample, pdf = sampler_dist_builder(x, y_corrected)

	k, mu, sigma = exponnorm.fit(sample)
	emg = exponnorm(k, mu, sigma) # initiate is only possible with parameters (different from sklearn) --> frozen distribution 
	lam = 1/(k*sigma)
	mean, var, skew, kurt = emg.stats(moments='mvsk') # provide moments (raw, central, standardized and standardized for m, v, s and k, respectively)
	# mean, var, skew, kurt = emg.moment(1), emg.moment(2), emg.moment(3), emg.moment(4)

	if plotting: 
		## sanity check
		# print(emg.moment(1)) # non-central 1st momemt
		# print(emg.moment(2)) # non-central 2nd momemt
		# print(emg.moment(3)) # non-central 3rd momemt 
		# print(emg.moment(4)) # non-central 4th momemt 
		# for param, name in zip([mu, sigma, k, lam], ["mu", "sigma", "k", "lam"]):print(f"{name:6} = {param:>8.5g}")# f-string padding '', '^' and '>' for left, center and right padding, number afterward determines padding size
		# print(mean)																	# 1st moment from stats() method
		# print(mu + (1/lam)) 															# 1st moment calculated from EMG parameters
		# print(var) 																	# 2nd moment from stats() method
		# print(sigma**2 + (1/(lam**2))) 												# 2nd moment calculated from EMG parameters
		# print(skew) 																	# 3rd moment from stats() method
		# print((2/((sigma**3) * (lam**3)))*((1+(1/((sigma**2) * (lam**2))))**(-3/2))) 	# 3rd moment calculated from EMG parameters
		# print((2/(sigma**3))) 	# 3rd moment calculated from EMG parameters

		## plotting for visual check
		##
		x_range = np.linspace(x.min(), x.max(), 100)
		plt.plot(x, pdf, "g:", lw=2, label="peak, scaled & \nbaseline-corrected,\nanalogous to 'true' PDF")
		plt.hist(sample, bins=60, density=True, alpha=0.4, label='resampled binned data')
		plt.plot(x_range, emg.pdf(x_range), 'r-', lw=2, alpha=0.6, label='EMG, fitted PDF')
		
		
		ax = plt.gca() # to access get_lines method
		line_count = len(ax.get_lines()) # gets the number of lines plotted with plt.plot(...
		if line_count < 1:
			color = "black"
		elif line_count > 2:
			color = "gray"
		elif line_count > 3:
			color = "silver"
		else:
			color = "black"
		plt.plot(x, (baseline - y + y_corrected), color = color, linestyle = (0, (5, 5)), label = "baseline")
		plt.legend(); plt.xlabel('x'); plt.ylabel('probability density');
		
	dist_moments_central: dict[str, float] = {"variance": var}
	dist_moments_standardized: dict[str, float] = {"skew": skew, "kurtosis": kurt}
	dist_moments_raw: dict [int, float] = {1: emg.moment(1), 2: emg.moment(2), 3: emg.moment(2), 4: emg.moment(4)}
	dist_params: dict[str, float] = {"k":k, "mu": mu, "sigma": sigma, "lambda": lam}

	return dist_moments_central, dist_moments_standardized, dist_moments_raw, dist_params, baseline


def _interp_crossing(x: np.ndarray, y: np.ndarray, level: int|float) -> tuple[float, float]:
	"""
	Find (interpolte) x-values corresponding to 'level' on both sides of the peak.
	Linear interpolation is used between the two closest 'y' points. If the peak data is noisy and crosses 'level' multiple times, the first pair of surronding points is selected for interpolation (first in the ascending half and first in the descending half). For very noisy data, smoothing before is recommened.

	Inputs:
		x (np.ndarray): x-coordinate data e.g. eluted volume
		y (np.ndarray): y-coordinate data e.g. Conductivity signal
		level (int|float): y-value to which interpolate the x-values

	Output:
		x_left, x_right: tu(list[float]|np.ndarray): baseline-corrected y-coordinate data e.g. baseline corrected signal
	"""

	idx_apex = np.argmax(y) # peak max

	## left half of the peak
	y_left = y[:idx_apex]
	idx0 = ((y_left - level)>=0).argmax()-1 # index "before" crossing 'level'
	idx1 = ((y_left - level)>=0).argmax() # index "after" crossing 'level'
	x_left = np.interp(level,y[idx0: idx1 + 1],x[idx0: idx1 + 1]) # interpolated value




	y_right = y[idx_apex:]
	idx0 = idx_apex + ((y_right - level)<=0).argmax()-1 # index "before" crossing 'level'
	idx1 = idx_apex + ((y_right - level)<=0).argmax() # index "after" crossing 'level'


	x_right = np.interp(level,y[idx0: idx1 + 1],x[idx0: idx1 + 1]) # interpolated value

	return x_left, x_right


def peak_asymmetry(x, y, percent_as=10, percent_tf=5) -> dict[str, int|float]:
    """
    Compute asymmmetry (asymmetry) at `percent_as`% and USP tailing factor (t_f) at `percent_tf`%.
    Requires y to be baseline-corrected
    
   Asymmetry factor at 10% height<br>
	$A_s = \frac{b}{a}$
	<br>with<br> 
	$a = x_{\text{apex}} - x_{L,0.10},\newline
	b = x_{R,0.10} - x_{\text{apex}}\newline$

	USP tailing factor at 5% height<br>
	$T_f = \frac{w_{0.05}}{2f}$
	<br>with<br> 
	$W_{0.05} = x_{R,0.05} - x_{L,0.05},\newline
	f = x_{\text{peak}} - x_{L,0.05}$

    
    Inputs: 
		x (np.ndarray): x-coordinate data e.g. eluted volume
        y (np.ndarray): y-coordinate data, baseline-corrected e.g. Conductivity signal
		baseline (None|int|float|np.ndarray): user-defined scalar, array (with same size as 'y') or min(y) if None
        percent_as = threshold level for Asymmetry calculation (default = 10%)
        percent_tf = threshold level for Tailing factor calculation (default = 5%)
        
    Outputs:
		asymmetry_dict (dict[str, float|int]): contains asymmetry and tailing factor data as well as intermediate calculations e.g. side a for asymmetry calculation
    """
    x = np.asarray(x); y = np.asarray(y)
    assert x.ndim == y.ndim == 1 and x.size == y.size and x.size >= 5
    
    # # create baseline array
    # if baseline is None:
    #     # Take the minimum y in the window as baseline (better: do a proper baseline correction upstream)
    #     correction_baseline = np.repeat(np.min(y), y.size)
    # elif np.isscalar(baseline):
    #     correction_baseline = np.repeat(float(baseline), y.size)
    # else:
    #     correction_baseline = np.asarray(baseline)
    #     if correction_baseline.shape != y.shape:
    #         raise ValueError("baseline array must match y length")
    
    # y_corr = y - correction_baseline

    # # peak apex
    # i_peak = np.argmax(y_corr)
    # x_peak, y_peak = x[i_peak], y[i_peak]
    # baseline_val = correction_baseline[i_peak]
    # h = y_peak
    # print(f"{y_peak=}")
    # print(f"{baseline_val=}")
    # print(f"{h=}")
    
    i_peak = np.argmax(y) # apex of peak, indice of
    x_peak, y_peak = x[i_peak], y[i_peak]
    h = y_peak
    


    if h <= 0:
        raise ValueError("peak height (h) is zero or lower, likely a problem with baseline")
    
    # target levels
    level_as = (percent_as/100.0) * h
    level_tf = (percent_tf/100.0) * h
    
    # ensure arrays are referenced to absolute signal (no need to shift, we just compare to levels)
    # find left/right intercepts by interpolation
    xL_as, xR_as = _interp_crossing(x, y, level_as)
    xL_tf, xR_tf = _interp_crossing(x, y, level_tf)
    
    # distances
    a = x_peak - xL_as
    b = xR_as - x_peak
    width_tf = xR_tf - xL_tf
    f = x_peak - xL_tf
    
    asymmetry = b / a if np.isfinite(a) and a > 0 else np.nan
    t_f = width_tf / (2*f) if np.isfinite(f) and f > 0 else np.nan
    
    asymmetry_dict= dict()
    asymmetry_dict.update({
        'As': asymmetry,
        'Tf': t_f,
        'xApex': x_peak,
        'yApex': y_peak,
        'A_10pct': a, 'B_10pct': b,
        'width_tf': width_tf, 'f': f})

    
    return asymmetry_dict


def cpt_emg_based(x_app: list[float]|np.ndarray, y_app: list[float]|np.ndarray, x_sys: list[float]|np.ndarray|None = None, y_sys: list[float]|np.ndarray|None = None, **kwargs):
	"""
	calculate column performance based on assumption that peak(s) follow Exponentially-modified Gaussian (EMG).
	Include option to correct for system contribution (i.e. extra-column contribution)

	Inputs:
		x (list[float]|np.ndarray): column peak x-coordinate data e.g. eluted volume
		y (list[float]|np.ndarray): column peak y-coordinate data e.g. Conductivity signal
		x_sys (list[float]|np.ndarray): system contribution x-coordinate data e.g. eluted volume
		y_sys (list[float]|np.ndarray): system contribution y-coordinate data e.g. Conductivity signal

		kwargs:
		plotting (bool): whether to plot (or not) the peak and EMG dist (PDF)
		baseline_method (str): baseline calculation method can be either 'spline' or 'interpolation'.
		col_bed_height (int|float): bed height in cm
		d_particle (int|float): resin particle diameter in µm
			
	Output:
		col_performance (dict[str, float]): dict with column performance data (hetp, reduced hetp, non-corrected skew) and column moments (corrected for extra-colum effect if provided)

	NOTES:
		system contribution refers to the injection without packed column to characterize the extra-column effects
	
	"""
	app_central_moments, app_std_moments, app_raw_moments, app_params, app_baseline = emg_fitter(x_app, y_app, **kwargs)
	# report apparent moments
	print("\tmoments (raw)")
	print("\t{0:>8} {1:>8} {2:>8} {3:>8} {4:>8}".format("", "1st", "2nd", "3rd", "4th"))
	print("\t{0:>8} {1:>8.4g} {2:>8.4g} {3:>8.4g} {4:>8.4g}".format("app", *app_raw_moments.values()))

	

	# ## draft
	# k, loc, scale, _ = apparent_moments.values()
	# # mom1, mom2, mom3 = 
	# print(f"\tnon-central 1st moment = {exponnorm(k,loc,scale).moment(1)}")
	# ## draft

	# if system contribution (extra-column) is available, correct for it
	if x_sys is not None and y_sys is not None:
		
		sys_central_moments, sys_std_moments, sys_raw_moments, sys_params, _  = emg_fitter(x_sys, y_sys, **kwargs)
		# report raw moments
		print("\t{0:>8} {1:>8.4g} {2:>8.4g} {3:>8.4g} {4:>8.4g}".format("sys", *sys_raw_moments.values()))

		try: # try moment-based devonvolution (moments are additive)
			col_central_moments = dict()
			for moment in [key for key in app_central_moments.keys() if key in sys_central_moments.keys()]:
				col_central_moments[moment] = app_central_moments[moment] - sys_central_moments[moment]
			col_raw_moments = dict()
			for moment in [key for key in app_raw_moments.keys() if key in sys_raw_moments.keys()]:
				col_raw_moments[moment] = app_raw_moments[moment] - sys_raw_moments[moment]
			col_std_moments = dict()			
			for moment in [key for key in app_std_moments.keys() if key in sys_std_moments.keys()]:
				col_std_moments[moment] = app_std_moments[moment] - sys_std_moments[moment]

			mean = col_raw_moments[1]
			var = col_central_moments["variance"]
			skew = col_std_moments["skew"]

			warnings.filterwarnings("error") # catch warning as exceptions
			lam = (2/skew)**(1/3)
			sigma = np.sqrt(var - (1/(lam**2)))
			mu = mean - (1/lam)
			k = 1/(sigma**lam)
			warnings.resetwarnings() # turn off warning-catching
			#print(mu, sigma, lam, k)

		except RuntimeWarning: 
			# negative 3rd moment results in imaginary lam, thus imossible. 
			# Some combinations of EMG distributions cannot be deconvoluted with momments 
			# (FYI convolution of two EMGs does not result in one EMG, but can in some cases be approcimated to)
			
			raise ValueError('negative 3rd moment results in imaginary lambda parameter for the distribution, thus imossible')
		
	else: 
		col_raw_moments.update(app_raw_moments)
		col_central_moments.update(app_central_moments)
		col_std_moments.update(app_std_moments)




	# report col raw moments
	print("\t{0:>8} {1:>8.4g} {2:>8.4g} {3:>8.4g} {4:>8.4g}".format("col", *col_raw_moments.values()))



	## column performance

	app_as = peak_asymmetry(x_app, y_app, baseline=app_baseline)

	col_bed_height = kwargs.get("col_bed_height", None)
	d_particle = kwargs.get("d_particle", None) # get/check if particle diameter is provided
 
	if col_bed_height is not None:
		var = col_central_moments["variance"]
		mean = col_raw_moments[1]
		hetp_cm = ((var * col_bed_height)/mean**2)
		hetp = hetp_cm * 1e+4
		print(f"\t{hetp = :.3g} \u03BCm") # "\u03BC" is for mu small greek letter

		if d_particle is not None:
			reduced_hetp = hetp / d_particle
			print(f"\t{reduced_hetp = :.3g}")

	print("\tAsymmetry = {:.3g} (non-corrected, apparent)".format(app_as["As"]))

	#return column performance
	col_performance = dict()
	col_performance.update({"hetp": hetp})
	col_performance.update({"reduced_hetp": reduced_hetp})
	col_performance.update({'skew_non_corrected': app_std_moments["skew"]})
	col_performance.update({"asymmetry" :app_as["As"]})
	col_performance.update(col_raw_moments)
	col_performance.update(col_central_moments)
	col_performance.update(col_std_moments)

	return col_performance


def moments_numerical_integration(x: np.ndarray, y: np.ndarray, k_max: int = 4, **kwargs) -> dict[str, dict]:
    """
    Compute statistical moments form 0 to k_max-th order for a single peak using numerical integration (trapezoidal rule).
	Parameters
    ----------
	x : np.ndarray
		1D array of x-coordinates (e.g., time).
	y : np.ndarray
		1D array of intensities/signal values (non-negative, baseline-corrected).
	n : int
		Maximum moment order to compute (>=2).
    kwargs : dict
		keyworded optional arguments
        
    Returns:
    ----------
    results_dict : dict[str, float]
     	dict containing the n-th order moments 
        'area': zeroth raw moment (∫ y dx)
    	'raw': normalized raw moments m_k = (∫ x^k y dx) / (∫ y dx)  for k=0..k_max
		'central': normalized central moments μ_k = (∫ (x-μ)^k y dx) / (∫ y dx) for k=0..k_max
		'standardized': skewness and kurtosis (excess)
      
    Notes:
      - Works for non-uniform x spacing.
      - Automatically sorts by x and clips tiny negative baselines.
      - Assumes a single, baseline-corrected peak in y ≥ 0.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], np.clip(y[order], 0.0, None)

    #plt.plot(x, y)

    # Zeroth moment raw moment, non-normalized (area)
    area = np.trapz(y, x)
    if area <= 0:
        raise ValueError("Area is zero or negative; check your data/baseline.")

    # Raw (about origin), normalized for peak area (total intensity)
    raw_non_norm = [np.trapz((x**k) * y, x) for k in range(k_max + 1)] # M_kmax
    raw_norm = [m_large / area for m_large in raw_non_norm]  # m0..m_kmax
    raw_norm = raw_non_norm / area # m0..m_kmax
    
	# Package results (in labeled dicts)
    mean = raw_norm[1] if k_max >= 1 else np.trapz(x * y, x) / area  # centroid
    raw_dict = {f"m_{k}": raw_norm[k] for k in range(k_max + 1)}

    # Central (about the mean)
    dx = x - mean # re-centering
    central_integrals = [np.trapz((dx**k) * y, x) for k in range(k_max + 1)]
    central = [central_int / area for central_int in central_integrals]  # μ0..μ_kmax
    # Package results in labeled dicts
    central_dict = {f"mu_{k}": central[k] for k in range(k_max + 1)}
    if central_dict["mu_1"] < 1e-9: # clip values approaching negative inf
        central_dict["mu_1"] = float(0) 
    
	# Standardized moments
    mu_2 = central_dict.get("mu_2", None)
    if mu_2 is None or mu_2 <= 0:
        raise ValueError("Variance (mu_2) must be positive to standardize moments.")
    standardized_dict = dict()
    for k in range(3, k_max + 1):  # standardized moments start at order 3
        gamma_k = central_dict[f"mu_{k}"] / (mu_2 ** (k / 2))
        standardized_dict[f"gamma_{k}"] = gamma_k

    results_dict = {
        "area": area,                       # ∫ y dx
        "raw": raw_dict,                    # normalized m0..m_kmax (m0 == 1)
        "central": central_dict,            # centered about the mean μ0..μ_kmax (μ0 == 1, μ1 == 0)
        "standardized": standardized_dict,	# standardized moments by the appropriate power of variance to standardized by disperison and be dimensionless 
    }

    return results_dict


def cpt_moment_numerical(
		x_app: list[float]|np.ndarray, 
		y_app: list[float]|np.ndarray, 
		x_sys: list[float]|np.ndarray|None = None, 
		y_sys: list[float]|np.ndarray|None = None, 
		**kwargs):
	"""
	Calculate column performance using statistical moments measured by numerical integration. 
	Includes option to correct for system contribution (i.e. extra-column contribution)

	Inputs:
		x (list[float]|np.ndarray): column peak x-coordinate data e.g. eluted volume
		y (list[float]|np.ndarray): column peak y-coordinate data e.g. Conductivity signal
		x_sys (list[float]|np.ndarray): system contribution x-coordinate data e.g. eluted volume
		y_sys (list[float]|np.ndarray): system contribution y-coordinate data e.g. Conductivity signal

		kwargs:
		plotting (bool): whether to plot (or not) the peak and EMG dist (PDF)
		baseline_method (str): baseline calculation method can be either 'spline' or 'interpolation'.
		col_bed_height (int|float): bed height in cm
		d_particle (int|float): resin particle diameter in µm
			
	Output:
		col_performance (dict[str, float]): dict with column performance data (hetp, reduced hetp, non-corrected skew) and column moments (corrected for extra-colum effect if provided)

	NOTES:
		system contribution refers to the injection without packed column to characterize the extra-column effects
	
	"""
	baseline_method: str = kwargs.get("baseline_method", "spline") # 'spline' as default methhod for baseline finding/correction
	y_app_corr, baseline_app = baseline_corrector(x_app, y_app, baseline_method = baseline_method)
	app_moments = moments_numerical_integration(x_app, y_app_corr, **kwargs)

	# if system contribution (extra-column) is available, correct for it
	col_moments = dict()
	if x_sys is not None and y_sys is not None:
		y_sys_corr, baseline_sys = baseline_corrector(x_sys, y_sys, baseline_method=baseline_method)
		sys_moments = moments_numerical_integration(x_sys, y_sys_corr, **kwargs)
		col_moments["raw"] = dict()
		col_moments["raw"]["m_0"] = 1
		col_moments["raw"]["m_1"] = app_moments["raw"]["m_1"] - sys_moments["raw"]["m_1"]
		col_moments["raw"]["m_2"] = app_moments["raw"]["m_2"] - sys_moments["raw"]["m_2"]
		col_moments["central"] = dict()
		col_moments["central"]["mu_2"] = app_moments["central"]["mu_2"] - sys_moments["central"]["mu_2"]
	else: 
		sys_moments = dict()
		col_moments["raw"] = dict()
		col_moments["raw"].update(app_moments["raw"])
		col_moments["central"] = dict()
		col_moments["central"].update(app_moments["central"])
		col_moments["standardized"] = dict()
		col_moments["standardized"].update(app_moments["standardized"])


	print("\traw moments (non-normalized):")
	print("\t{0:>8} {1:>8}{2}".format("", "0th", " (peak area)")) 
	print("\t{0:>8} {1:>8.4g}".format("app", app_moments["area"]))
	try: 
		print("\t{0:>8} {1:>8.4g}".format("sys", sys_moments["area"]))
	except KeyError:
		pass

	print("\traw moments (normalized)")
	print("\t{0:>8} {1:>8} {2:>8} {3:>8} {4:>8} {5:>8}".format("", "0th", "1st", "2nd", "3rd", "4th")) 
	try: # case were x_sys and y_sys (extra-column effects, system contributions) are not provided
		print("\t{0:>8} {1:>8.4g} {2:>8.4g} {3:>8.4g} {4:>8.4g} {5:>8.4g}".format("col", *col_moments["raw"].values(), "-" ))
	except ValueError: # case where the col_moments are calculated removing the system contribution
		print("\t{0:>8} {1:>8.4g} {2:>8.4g} {3:>8.4g} {4:>8} {5:>8}".format("col", *col_moments["raw"].values(), "-", "-" ))
	print("\t{0:>8} {1:>8.4g} {2:>8.4g} {3:>8.4g} {4:>8.4g} {5:>8.4g}".format("app", *app_moments["raw"].values()))
	try: 
		print("\t{0:>8} {1:>8.4g} {2:>8.4g} {3:>8.4g} {4:>8.4g} {5:>8.4g}".format("sys", *sys_moments["raw"].values()))
	except KeyError:
		pass

	## column performance
	as_app = peak_asymmetry(x_app, y_app_corr, baseline=baseline_app)

	col_bed_height = kwargs.get("col_bed_height", None)
	d_particle = kwargs.get("d_particle", None) # get/check if particle diameter is provided

	if col_bed_height is not None:
		var = col_moments["central"]["mu_2"]
		mean = col_moments["raw"]["m_1"]
		hetp_cm = ((var * col_bed_height)/mean**2)
		hetp = hetp_cm * 1e+4
		print(f"\t{hetp = :.3g} \u03BCm") # "\u03BC" is for mu small greek letter

		if d_particle is not None:
			reduced_hetp = hetp / d_particle
			print(f"\t{reduced_hetp = :.3g}")

	print("\tAsymmetry = {:.3g} (non-corrected, apparent)".format(as_app["As"]))
	print("\tskew = {:.3g} (non-corrected, apparent)".format(app_moments["standardized"]["gamma_3"]))

	## return column performance metrics	
	col_performance = dict()
	col_performance["performance"] = dict()
	col_performance["performance"].update({"hetp": hetp})
	col_performance["performance"].update({"reduced_hetp": reduced_hetp})
	col_performance["performance"].update({'skew_non_corrected': app_moments["standardized"]["gamma_3"]})
	col_performance["performance"].update({"asymmetry_non_corrected" :as_app["As"]})
	
	## add the remaining dicts with moments (apparent, column, system)
	col_performance["app_moments"] = dict()
	col_performance["app_moments"].update(app_moments)
	col_performance["col_moments"] = dict()
	col_performance["col_moments"].update(col_moments)
	col_performance["sys_moments"] = dict()
	col_performance["sys_moments"].update(sys_moments)
		
	return col_performance


# ─────────────────────────────────────────────────────────────────────────────
# 1.  FUNCTION LIBRARY (with docstrings & type annotations)
# ─────────────────────────────────────────────────────────────────────────────
def savgol_smooth(
    y: np.ndarray,
    *,
    points: int | None = None,
    width: float | None = None,
    dx: float = 1.0,
    polyorder: int = 3,
    deriv: int = 0,
    mode: str = "interp") -> np.ndarray:
	
    """
    Savitzky–Golay *convolution-based* smoother with two
    alternative ways to choose the window size.

    Parameters
    ----------
    y : ndarray
        1-D signal to be smoothed.
    points : int, optional
        Explicit odd window length *in samples*.
        Exactly one of ``points`` or ``width`` must be given.
    width : float, optional
        Desired physical window width (same units as ``dx``); the
        function converts it to an odd number of samples automatically.
    dx : float, optional
        Sample spacing used when ``width`` is supplied.  Default is 1.0.
    polyorder : int, optional
        Polynomial order for the local least-squares fit
        (default *3*).  Must satisfy ``polyorder < window_length``.
    deriv : int, optional
        Order of the derivative to return  
        (0 ⇒ smoothed signal, 1 ⇒ first derivative, …).
    mode : str, optional
        End-point handling; forwarded to
        :func:`scipy.signal.savgol_filter`.  Default ``"interp"`` gives
        the smoothest edge behaviour.

    Returns
    -------
    y_smooth : ndarray
        Array of the same shape as `y` containing either the smoothed
        signal (`deriv = 0`) or its requested derivative.

    Notes
    -----
    *   Set either ``points`` **or** ``width``.  If both (or neither) are
        supplied a ``ValueError`` is raised.
    *   For chromatograms sampled at irregular *x* positions, supply a
        representative ``dx`` (e.g. the median spacing) when using
        ``width``.
    *   The derivative option lets you do peak-finding or baseline
        detection without re-computing the filter coefficients.
    """
    # ── argument validation ────────────────────────────────────────────────
    if (points is None) == (width is None):
        raise ValueError("Specify exactly one of `points` or `width`.")

    if width is not None:                       # convert width → points
        points = int(np.floor(width / dx)) | 1  # enforce odd using bit-or 1
        if points < polyorder + 2:              # widen if too small
            points = polyorder + 2 | 1

    if points % 2 == 0:
        points += 1  # force odd

    if polyorder >= points:
        raise ValueError("`polyorder` must be less than `window_length`.")

    # ── apply Savitzky–Golay filter ───────────────────────────────────────
    y_smooth: np.ndarray = savgol_filter(
        y, window_length=points, polyorder=polyorder,
        deriv=deriv, delta=dx, mode=mode
    )
    return y_smooth


def quantile_bounds(x: np.ndarray, y: np.ndarray, tol: float = 1e-8) -> tuple[float, float]:
    """
    Return lower/upper *x* positions that enclose (1 - 2 x tol) of the peak area.

    Parameters
    ----------
    x, y : ndarray
        Abscissae and positive ordinates of a 1-D peak.
    tol : float, optional
        Tail probability cut-off on each side (default ``1e-8``).

    Returns
    -------
    x_low, x_high : float
        *x* coordinates where the cumulative area reaches
        ``tol`` and ``1 - tol`` respectively.
    """
    cdf: np.ndarray = integrate.cumulative_trapezoid(y, x, initial=0)
    cdf /= cdf[-1]
    x_low = float(np.interp(tol, cdf, x))
    x_high = float(np.interp(1 - tol, cdf, x))
    return x_low, x_high


def build_grid(
    data: list[tuple[np.ndarray, np.ndarray]],
    n: int = 2**14,
    tol: float = 1e-8,
    pad_factor: int = 4,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Construct an equally spaced *core* grid that covers **all** input peaks
    up to the requested quantile tolerance and append zero-padding.

    Parameters
    ----------
    data : list of (x, y) tuples
        Each tuple contains abscissæ and ordinates of a curve.
    n : int, optional
        Number of points in the core grid (power of two recommended).
    tol : float, optional
        Tail probability passed to :func:`quantile_bounds`.
    pad_factor : int, optional
        Overall array length multiplier after padding
        (``pad_factor = 4`` ⇒ 3 × n zeros added).

    Returns
    -------
    x_core : ndarray, shape (n,)
        Evenly spaced grid spanning the informative support.
    x_full : ndarray, shape (pad_factor · n,)
        ``x_core`` followed by zeros — ready for FFT algorithms that assume
        periodicity.
    dx : float
        Grid spacing, useful for normalisation or integration.
    """
    lows, highs = zip(*(quantile_bounds(x, y, tol) for x, y in data))
    x_core: np.ndarray = np.linspace(min(lows), max(highs), n)
    x_full: np.ndarray = np.concatenate([x_core, np.zeros((pad_factor - 1) * n)])
    dx: float = float(x_core[1] - x_core[0])
    return x_core, x_full, dx


def interp_to_grid(
    x_src: np.ndarray, y_src: np.ndarray, x_target: np.ndarray
) -> np.ndarray:
    """
    Linearly interpolate ``(x_src, y_src)`` onto ``x_target``.

    Values outside the original ``x_src`` range are set to zero.

    Returns
    -------
    y_interp : ndarray, shape = ``x_target.shape``
        Interpolated ordinates on the new grid.
    """
    f = interpolate.interp1d(
        x_src, y_src, kind="linear", bounds_error=False, fill_value=0.0
    )
    return f(x_target)


def richardson_lucy(
    obs: np.ndarray, kernel: np.ndarray, iters: int = 300
) -> np.ndarray:
    """
    1-D Richardson-Lucy deconvolution (non-negativity-preserving ).

    Parameters
    ----------
    obs : ndarray
        Observed signal *A* (must be non-negative) on an equi-spaced grid.
    kernel : ndarray
        Convolution kernel *S* on the same grid (will be internally
        normalised to unit sum).
    iters : int, optional
        Number of RL iterations (default 300).  
        More iterations sharpen the result but may amplify noise.

    Returns
    -------
    est : ndarray
        Deconvolved estimate of the hidden signal *C* on the same grid.
    """
    kernel = kernel / kernel.sum()
    kernel_flip = kernel[::-1]
    est = np.full_like(obs, obs.mean())  # positive initial guess

    for _ in range(iters):
        conv = fft.ifft(fft.fft(est) * fft.fft(kernel)).real
        conv[conv <= 0] = 1e-12  # avoid division by zero
        ratio = obs / conv
        est *= fft.ifft(fft.fft(ratio) * fft.fft(kernel_flip)).real
        est[est < 0] = 0  # enforce non-negativity

    return est
