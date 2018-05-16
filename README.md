# Updated_T-LoCoH_Algorithm

The algorithm, as presented in Movement Ecology (https://doi.org/10.1186/s40462-017-0110-4), carries out a grid-based search of k and s values in an attempt to optimize the T-LoCoH parameter selection process. This particular approach values consistency in home range delineations (as measured through cross-validation of alternative training-testing sets) over alternatives, such as avoiding a swiss-cheese appearance.

The T-LoCoH_Algorithm_Functions_Efficient.R file consists of a set of functions that are explained in detail here:

1) format.data - This function takes a data frame ('data') and a character string ('crs') that identifies the target projection. The input data frame is required to consist of at least three columns called 'long' (for longitude), 'lat' (for latitude), and 'Datetime' (character string or POSIXct of the form yyyy-mm-dd hh:mm:ss). The output is a dataframe with corresponding 'x', 'y', and 'datetime' columns, but with any NA points filled using a Kalman Smoother.

2) train.test - This function takes the output data frame from format.data as well as an optional seed value for the sake of reproducibility (set at 1 as the default). The result is a new data frame with 100 columns and rows equal to the length of the movement track. Each column represents a different training-testing split consisting of approximately 0.22% testing points, approximately 20% removed to reduce the effect of autocorrelation, and approximately 80% training data upon which the hullsets are constructed.

3) find.test.pts - This function takes the output of train.test and projects the test points as SpatialPoints. It is called automatically within the next function ('algo.efficient')

4) algo.efficient - This function takes five arguments:
 - a) tt.split - the result of the train.test function
 - b) k.max - the maximum value of k that the search will reach
 - c) data.lxy - a properly formatted lxy-object from tlocoh
 - d) data - the original dataset (xy-coords are sufficient)
 - e) crs - the a character string of the projection in proj4string form

This command uses a parallelized loop to perform a set of three grid based searches that increase in resolution. The first searches between 20 and k.max in intervals of 20 (but also including 4 as the smallest k value) and between 0 and 0.05 in increments of 0.1. The second grid search takes the optimal k value (opt1) from the previous loop and searches over the k values from opt1 - 20 to opt1 + 20 in increments of 5 and between 0 and 0.05 in increments of 0.1. Finally, the optimal k value from that search is identified (opt2) and an even smaller range (from opt2 - 5 to opt2 + 5) is searched in increments of 1 and over s values from 0 to 0.05 in increments of 0.001. The entire set of values emerging from this search is returned and the optimal parameter set can be identified from that.

The code here is also implemented by looping over a set of paths (name.list) arbitrariliy named c('AA', 'BB', 'CC', 'DD'). The user can simply replace these with their .csv file names (excluding the '.csv' portion) to run the efficient algorithm on their data. The result of the loop implemented here is a single data frame (opt.params) that has a single row for each individual with the optimal parameter values, the total area of the optimal hullset, and the cumulative probability as calculated following the equation in the corrected article.
