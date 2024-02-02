#!/usr/bin/env python
""" 
Python script to download selected files from rda.ucar.edu.
After you save the file, don't forget to make it executable
i.e. - "chmod 755 <name_of_script>"
"""
import sys, os
from urllib.request import build_opener

opener = build_opener()

filelist = [
  'https://stratus.rda.ucar.edu/ds083.2/grib2/2023/2023.08/fnl_20230815_00_00.grib2',
  'https://stratus.rda.ucar.edu/ds083.2/grib2/2023/2023.08/fnl_20230815_06_00.grib2',
  'https://stratus.rda.ucar.edu/ds083.2/grib2/2023/2023.08/fnl_20230815_12_00.grib2',
  'https://stratus.rda.ucar.edu/ds083.2/grib2/2023/2023.08/fnl_20230815_18_00.grib2',
  'https://stratus.rda.ucar.edu/ds083.2/grib2/2023/2023.08/fnl_20230816_00_00.grib2',
  'https://stratus.rda.ucar.edu/ds083.2/grib2/2023/2023.08/fnl_20230816_06_00.grib2',
  'https://stratus.rda.ucar.edu/ds083.2/grib2/2023/2023.08/fnl_20230816_12_00.grib2',
  'https://stratus.rda.ucar.edu/ds083.2/grib2/2023/2023.08/fnl_20230816_18_00.grib2'
]

for file in filelist:
    ofile = os.path.basename(file)
    sys.stdout.write("downloading " + ofile + " ... ")
    sys.stdout.flush()
    infile = opener.open(file)
    outfile = open(ofile, "wb")
    outfile.write(infile.read())
    outfile.close()
    sys.stdout.write("done\n")
