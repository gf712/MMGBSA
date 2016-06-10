## MMGBSA parser for AMBER 15

* ### mmgbsa_parser.sh
Shell script to extract data from AMBER 15 MMGBSA.py output files.
Add this script to your path, so that it can be called by the python code. <br /> 
Future plans include writing the code in perl to speed up the data extraction.

* ### MMGBSA_parser.py
Python code to be excuted from the command line.
For example: <br /> 

```
MMGBSA_parser.py -i ./input_dir -o ./output_dir -fo test_file -pt MMGBSA test
```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Leading to the following plot in ```./output_dir```: <br /> 
![alt tag](https://raw.githubusercontent.com/gf712/MMGBSA/master/MMGBSA%20test.png)