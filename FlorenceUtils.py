
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

	#print sum(red_log_vector)

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

def load_arm1_data(file_name):

    curve_data = pd.read_pickle(file_name)
    arm1_data = curve_data[0:4]

    return arm1_data

def exp_func(x, a, b):

    return a*np.exp(-b*x)

def find_largest_error(data, organ, log_array):
    
    """Function to return the largest array for a given organ for a 
    given set of removed time points given by the @log_array"""
    
    error_vector = np.ones(4)

    
    if organ == 'WB':
        
        do_wb = True
        
    else:
        
        do_wb = False
    
    for i in range(4):
        
        t_vec = data['Time'][i]
        a_vec = data[organ][i]
        
        error_vector[i] = FlorenceUtils.return_diff_full_and_reduced(t_vec, a_vec, log_array, is_wb=do_wb)

    #print error_vector
    larges_error = max(np.min(error_vector), np.max(error_vector), key = abs)
    
    return larges_error
    

def construct_error_matrix(data, cut_matrix):
    
    """
    Should return the largest error for six different cut vectors for all
    pre-defined organs
    
    """
    
    error_matrix = np.zeros([6,4])
    
    organ_names = ['Liver', 'Spleen', 'Kidney', 'WB']
    
    for j in range(4):
    
        for i in range(6):
        
            cut_vec = cut_matrix[i,:]
            largest_organ_error = find_largest_error(data, organ_names[j], cut_vec)
        
            error_matrix[i,j] = largest_organ_error


    return error_matrix*100

def on_click(event):
    
    #print "pressed"
    
    #print len(polygons)
        
    if polygon_draw.contains_point((event.x, event.y)):
        
        cut_matrix = logical_back
        error_matrix = construct_error_matrix(data, cut_matrix)
        render_heatmap(error_matrix, cut_matrix)
        
    if polygon_restart.contains_point((event.x, event.y)):
        
        print "Restarting..."
        plt.clf()
        
    
    for i in range(0, 36): # TODO: Remove magic number...
                
        curr_pol = polygons[i]
                
        if curr_pol.contains_point((event.x, event.y)):
            
            curr_pol.set_facecolor('#404040')
                
            curr_pol_ind = poly_inds[i]
            
            logical_back[curr_pol_ind[1], curr_pol_ind[0]] = 1
            
                
    fig.canvas.draw()

def render_heatmap(error_matrix, cut_matrix = 0):
    
    """Function to render the heat map from the error matrix"""
    
    copy = error_matrix
    chararray = error_matrix_to_annot(error_matrix)
    y_label = cut_points_string_from_matrix(cut_matrix)
    x_label = ['Liver', 'Spleen', 'Kidney', 'WB']
    
    fig = plt.figure()
    g = sns.heatmap(error_matrix, 
                annot=chararray,
                fmt='', 
                cbar = False, 
                square=True,
                vmax = 10,
                vmin = -10,
                yticklabels = y_label,
                xticklabels = x_label,
                )
        
    g.set_yticklabels(g.get_yticklabels(), rotation = 0)
    
    return fig

def error_matrix_to_annot(error_matrix):
    
    """Utility function to convert the error function into
       a matrix of strings that can be read by sns.heatmap"""
    
    error_matrix[error_matrix > 10] = 10
    error_matrix[error_matrix < -10] = 10

    # Start of magic
    string_arr = pd.DataFrame(error_matrix/100).applymap(lambda x: '{:.1%}'.format(x)).values
    string_arr[string_arr == '10.0%'] = '>10%'
    # End of magic DO NOT TOUCH!
    
    return string_arr

def cut_points_string_from_matrix(cut_matrix):
    
    """ Utility function that handles the legend of the heatmap
    
        First assume that we have six different time points containing
        six possible time points each
        
    """
    
    glob_time = np.array(['2,', '4,', '8,', '24,', '96,', '168 '])
    
    time_cuts = []
    
    for i in range(6):
        
        cut_line = cut_matrix[i,:]
        I = [cut_line == 1]
        
        string_time = str(glob_time[I])
        
        string_time = string_time.replace("'", "")
        #string_time = string_time.replace("[", "")
        #string_time = string_time.replace("]", "")
        
        time_cuts.append(string_time)
        
        
    return time_cuts
