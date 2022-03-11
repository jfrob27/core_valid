# LHYRICA's core validation

This application allows you to validate cores catalogues in order to train LHYRICA's neural network.

## Installation

```
git clone https://github.com/jfrob27/core_valid.git
```

## Requirements

First, we recommend to create a python virtual environment:

```
cd core_valid
python3 -m venv lhyrica-env
source lhyrica-env/bin/activate
pip install -r requirements.txt
```

##Download the data

The next step is to download the data, the fits maps and the core catalogue to validate. 

## Launch Bokeh server

The core validation interface is based on [Bokeh](https://docs.bokeh.org/en/latest/index.html) Python library for interactive visualisations. Once all the requirements are installed, in the terminal in the *core_valid* folder, with the virtual environment activated, enter the following command:

```
bokeh serve --show core_validation.py
```

Then a window will automatically open in your web browser. Wait a few seconds to in order to charge the image and the current catalogue.