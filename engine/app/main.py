import os
import sys
import uuid
import json
import numpy as np
import pandas as pd
from flask import Flask, session, request, redirect, url_for, jsonify, render_template, send_from_directory
from sklearn.manifold import TSNE

APP = Flask(__name__)

@APP.route('/', methods=['GET', 'POST'])
def main_page():

	#
	# Data
	#

	with open('static/results/features.json', 'r') as fj:
		DATA = json.load(fj)

	Dict = {'data': DATA['DATA'], 'labx': DATA['LABS'][0], 'laby': DATA['LABS'][1],'name': {'name': 'QCB Analysis'}}

	return render_template('index.html', PythonVars=Dict)

if __name__ == '__main__':

	APP.secret_key = 'A0Zr98j/3yX R~XHH!jmN]LWX/,?RT'

	#APP.run(host='10.128.62.73', port=8080, debug=True)
	APP.run(port=8080, debug=True)
 


