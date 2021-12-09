 !# Multiples: 0 for none, 1 for Moho, 2 for all first-order
	integer, parameter :: mults = 2  
 !# Number of samples per trace
	integer, parameter :: nsamp = 512
 !# Sample rate (seconds)
	real, parameter :: dt = 0.033 
 !# Gaussian pulse width (seconds)! a =sqrt(2)/w  or w =sqrt(2)/a
	real, parameter :: width = 0.1 
 !# Alignment: 0 is none, 1 aligns on primary phase (P or S)
	integer, parameter :: align =1
 !# Shift of traces -- t=0 at this time (sec)
	real, parameter :: shift = 1
 !# Rotation to output: 0 is NS/EW/Z, 1 is R/T/Z, 2 is P/SV/SH
	integer, parameter :: out_rot = 1 
 !# Noise added to data (1 is max amplitude of the p-wave)
	real, parameter :: noise = 0