# Example usage
The json files are in Stats folder. The code below demonstrate how to draw plots with the json files.
``` sh
git clone git@github.com:Hsin-Yeh/SNSPD_Analysis.git
cd SNSPD_Laser_measurements

# Example usage
# Sweep photon number
python3 python/NIPXIE/plot_examples.py Stats/4p68K/*_1000mV_240degrees_4p68K_stats.json
python3 python/NIPXIE/plot_examples.py Stats/11p47K/*_104mV_240degrees_11p47KK_stats.json
# Sweep bias current
python3 python/NIPXIE/plot_examples.py Stats/4p68K/1011_*mV_240degrees_4p68K_stats.json
```
# Keys stored in the stat dict: 

``` python
# Analysis parameter
'pulse_integral', 'eff', 'pulse_range', 'pulse_range_err', 'pulse_range_ptp', 'pulse_range_ptp_err', 'pre_range', 'pre_range_err', 'pulse_fall_tau', 'avgMax', 'avgMin', 'avg', 'range_avg', 'resist'
# Probably the most useful will be 'pulse_range_ptp' and 'pulse_fall_tau'

# Spectrums
'spectrum'

# Histograms
'h_pulse_fall_range_ptp', 'h_pulse_fall_tau', 'h_pulse_FWHM'])
```


# Figures

They are stored in the plots folder
