#!/bin/sh

cp cme.ini cme.ini.OLD

rm -f out/test_suite/minimal_wave_angle_230/*
cp in/test_suite/minimal_wave_angle_230/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_wave_angle_270/*
cp in/test_suite/minimal_wave_angle_270/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_wave_angle_310/*
cp in/test_suite/minimal_wave_angle_310/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_with_intervention_wave_angle_215/*
cp in/test_suite/minimal_with_intervention_wave_angle_215/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_with_intervention_wave_angle_235/*
cp in/test_suite/minimal_with_intervention_wave_angle_235/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_with_intervention_wave_angle_270/*
cp in/test_suite/minimal_with_intervention_wave_angle_270/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_with_intervention_wave_angle_285/*
cp in/test_suite/minimal_with_intervention_wave_angle_285/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_with_intervention_wave_angle_305/*
cp in/test_suite/minimal_with_intervention_wave_angle_305/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_with_two_wavestations/*
cp in/test_suite/minimal_with_two_wavestations/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_with_sediment_input/*
cp in/test_suite/minimal_with_sediment_input/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/Happisburgh/*
cp in/test_suite/Happisburgh/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/Manuel_C003_0001/*
cp in/test_suite/Manuel_C003_0001/cme.ini .
./cme
echo ===============================================================================

mv cme.ini.OLD cme.ini
