#!/bin/sh

cp cme.ini cme.ini.OLD

# rm -f out/test_suite/GMD2017/CliffFineBays/*
# cp in/test_suite/GMD2017/CliffFineBays/cme.ini .
# ./cme
# echo ===============================================================================

rm -f out/test_suite/GMD2017/Groyne/*
cp in/test_suite/GMD2017/Groyne/cme.ini .
./cme
echo ===============================================================================

# rm -f out/test_suite/GMD2017/CliffFineSandBays/*
# cp in/test_suite/GMD2017/CliffFineSandBays/cme.ini .
# ./cme
# echo ===============================================================================

rm -f out/test_suite/GMD2017/UndefendedCoastline/*
cp in/test_suite/GMD2017/UndefendedCoastline/cme.ini .
./cme
echo ===============================================================================

mv cme.ini.OLD cme.ini
