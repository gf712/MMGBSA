## MMGBSA parser for AMBER 15

Python code to be executed from the command line.

Usage info:
```
$ MMGBSA_parser.py  --help
usage: /Users/je714/Projects/MMGBSA/MMGBSA_parser.py

optional arguments:
  -h, --help            show this help message and exit
  -inf INFO_FILE, --info_file INFO_FILE
                        The information file that is printed after the MMPBSA
                        calculation. Default is _MMPBSA_info
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory for all the generated files.
  -fo OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file name.
  -pt PLOT_TITLE, --plot_title PLOT_TITLE
                        Plot title.
  -ts TIME_STEP, --time_step TIME_STEP
                        Time step (in ns) between frames
  -v VERBOSE, --verbose VERBOSE
                        Switch verbose on/off. Default is True.
  -p PLOT, --plot PLOT  Switch plotting on/off. Default is True.

Script to extract information from MMGBSA.py Amber 15 output.
```
For example:

```bash
MMGBSA_parser.py -inf _MMPBSA_info
```

Leading to the following plot in ```./output_dir/plot.pdf```:

![alt tag](http://imgur.com/xKJhxSC)
