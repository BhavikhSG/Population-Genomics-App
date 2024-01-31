from flask import Flask, render_template
import pandas as pd
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
app = Flask(__name__)

app.config['SECRET_KEY'] = '78y23hryufsbu!kjnf7&'
@app.route('/')
def index():
    return render_template("index.html")

@app.route('/AnalysisType')
def AnalysisType():
    return render_template("Analysis.html")

if __name__ == '__main__':
    app.run(debug=True)
