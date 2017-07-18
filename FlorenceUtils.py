
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

def integrate_box_and_tail(t_series, a_series):

    # First the initial box

    box_area = t_series[0]*a_series[0]

    # Now the trapezoid

    I = 0

    for i in range(0, len(a_series)-2):

	I = I + ((a_series[i] + a_series[i + 1])/2)*(t_series[i+1]-t_series[i])

    A = box_area + I
    
    # Last the last two points, use exponential 'fit' from two points

    lam = (np.log(a_series[-1])-np.log(a_series[-2]))/(t_series[-1]-t_series[-2])

    T = a_series[-2]/lam*-1

    A = A + T
    
    return A

def return_diff_full_and_reduced(t_series, a_series, red_log_vector, is_wb = False):

    I_full = integrate_box_and_tail(t_series, a_series)

    index_delete = np.where(red_log_vector == 0)[0]

    t_red = np.delete(t_series, index_delete)
    a_red = np.delete(a_series, index_delete)

    if len(a_red) < 3 or len(t_red) < 3:

	print "Warning - series is too short, default to none"
	return None

    if is_wb:

	print sum(red_log_vector)

	if sum(red_log_vector) == 2 or sum(red_log_vector) == 3:

	    initial = np.array([1000,-0.001])
	    popt_red, pcov_red = curve_fit(exp_func, t_red, a_red, p0 = initial)

	else:
	    
	    popt_red, pcov_red = curve_fit(exp_func, t_red, a_red)

	popt, pcov = curve_fit(exp_func, t_series, a_series)
	
	I_full = a_series[0]/popt[1]
	I_red = a_series[0]/popt_red[1]

	diff = (I_red-I_full)/I_full

    else:

	I_red = integrate_box_and_tail(t_red, a_red)

    	diff = (I_red-I_full)/I_full

    return diff

def return_cut_string(glob_string, red_log_vector):

    index_delete = np.where(red_log_vector == 0)[0]

    glob_string_red = np.delete(glob_string, index_delete)

    return glob_string_red

def return_error_array(data, organ_list, logical_array):

    """Legacy function, replace by better function"""

    num_of_cuts = 6
    
    largest_errors = np.zeros([num_of_cuts, len(organ_list)])
    
    for k in range(0, num_of_cuts):
    
        largest_val_list = np.zeros(len(organ_list))
        
        for j in range(0, len(organ_list)):
        
            organ = organ_list[j]
            largest_val_organ = -1
        
            for i in range(0, 4):
            
                a_series = data[organ][i]
                t_series = data['Time'][i]

		_error =  return_diff_full_and_reduced(t_series, a_series, logical_array[k,:])

		if _error == None:

		    print "Invalid value encountered"

		else:

		    percentage_error = np.abs(_error)*100
        
		    if percentage_error > largest_val_organ:
        
        	        largest_val_organ = percentage_error
        
            largest_val_list[j] = largest_val_organ
    
        largest_errors[k,:] = largest_val_list

    return largest_errors



def load_arm1_data(file_name):

    curve_data = pd.read_pickle(file_name)
    arm1_data = curve_data[0:4]

    return arm1_data

def exp_func(x, a, b):

    return a*np.exp(-b*x)
