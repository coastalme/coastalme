set datafile separator ','

set key autotitle columnhead
set ylabel "Z (m)"
set xlabel 'Distance (m)'

plot "../out/test_suite/minimal_cons_wave_angle_230/profile_015_timestep_1200.csv" using 1:4 with lines
