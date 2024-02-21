from flask import Flask, render_template, request, jsonify
import pandas as pd
import os
from flask_sqlalchemy import SQLAlchemy
import plotly.graph_objs as go  # Import Plotly library
import sqlite3
import numpy as np
from sklearn.decomposition import PCA
from sqlalchemy import func



app = Flask(__name__)
app.config['SECRET_KEY'] = '78y23hryufsbu!kjnf7&'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///C:/Users/xavia/OneDrive/Documents/Masters-PG/Bioinformatics MSc/Software Development Group Project/WebApp/instance/DatabaseTest 3.db' ##configures the URI for connecting to the db; specifies path to the SQLite db file
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SQLALCHEMY_ECHO'] = True



db = SQLAlchemy(app) ##Initialise SQLAlchemy; creates SQLAlchemy object 'db' that is bound to the Flask app; used to define db models + perform db operations


##Database migration: Database models; Defines several classes that inherit from 'db.Model'=base class provided by SQLAlchemy for defining db models. each class=table in db. each attribute=column in table.
class PCA(db.Model):
    __tablename__ = 'PCA'
    SAMPLE_ID = db.Column(db.Integer, primary_key=True)
    PC1 = db.Column(db.Float)
    PC2 = db.Column(db.Float)
    POPULATION = db.Column(db.String, db.ForeignKey('POP_GROUP.POPULATION_ID'))
    SUPERPOPULATION = db.Column(db.String, db.ForeignKey('SUPERPOP.SUPERPOPULATION_ID'))

class POP_GROUP(db.Model):
    __tablename__ = 'POP_GROUP'
    POPULATION_ID = db.Column(db.String, primary_key=True, unique=True)
    POPULATION_NAME = db.Column(db.String)
    pca = db.relationship('PCA', backref='population', lazy=True)

class SUPERPOP(db.Model):
    __tablename__ = 'SUPERPOP'
    SUPERPOPULATION_ID = db.Column(db.String, primary_key=True, unique=True)
    SUPERPOPULATION_NAME = db.Column(db.String)

class VARIANTS(db.Model):
    __tablename__ = 'VARIANTS'
    CHROMO = db.Column(db.String)
    SNP_ID = db.Column(db.String, primary_key=True)
    POS = db.Column(db.Integer)
    REF = db.Column(db.String)
    ALT = db.Column(db.String)
    GENE = db.Column(db.String)
    CLINICAL_RELEVANCE = db.Column(db.String)

class AlleleFrequency(db.Model):
    __tablename__ = 'ALLELE_FREQUENCY'
    SNP_ID = db.Column(db.String, db.ForeignKey('VARIANTS.SNP_ID'), primary_key=True)
    ACB = db.Column(db.String)
    ASW = db.Column(db.String)
    BEB = db.Column(db.String)
    CDX = db.Column(db.String)
    CEU = db.Column(db.String)
    CHB = db.Column(db.String)
    CHS = db.Column(db.String)
    CLM = db.Column(db.String)
    ESN = db.Column(db.String)
    FIN = db.Column(db.String)
    GBR = db.Column(db.String)
    GIH = db.Column(db.String)
    GWD = db.Column(db.String)
    IBS = db.Column(db.String)
    ITU = db.Column(db.String)
    JPT = db.Column(db.String)
    KHV = db.Column(db.String)
    LWK = db.Column(db.String)
    MSL = db.Column(db.String)
    MXL = db.Column(db.String)
    PEL = db.Column(db.String)
    PJL = db.Column(db.String)
    PUR = db.Column(db.String)
    SIB = db.Column(db.String)
    STU = db.Column(db.String)
    TSI = db.Column(db.String)
    YRI = db.Column(db.String)

class GenotypeFrequencies(db.Model):
    __tablename__ = 'GENOTYPE_FREQUENCIES'
    SNP_ID = db.Column(db.String, db.ForeignKey('VARIANTS.SNP_ID'), primary_key=True)
    ACB = db.Column(db.String)
    ASW = db.Column(db.String)
    BEB = db.Column(db.String)
    CDX = db.Column(db.String)
    CEU = db.Column(db.String)
    CHB = db.Column(db.String)
    CHS = db.Column(db.String)
    CLM = db.Column(db.String)
    ESN = db.Column(db.String)
    FIN = db.Column(db.String)
    GBR = db.Column(db.String)
    GIH = db.Column(db.String)
    GWD = db.Column(db.String)
    IBS = db.Column(db.String)
    ITU = db.Column(db.String)
    JPT = db.Column(db.String)
    KHV = db.Column(db.String)
    LWK = db.Column(db.String)
    MSL = db.Column(db.String)
    MXL = db.Column(db.String)
    PEL = db.Column(db.String)
    PJL = db.Column(db.String)
    PUR = db.Column(db.String)
    SIB = db.Column(db.String)
    STU = db.Column(db.String)
    TSI = db.Column(db.String)
    YRI = db.Column(db.String)


##Flask routes to handle HTTP requests:

@app.route('/') ##'/' route = index pg
def index():
    populations = db.session.query(PCA.POPULATION).distinct().all() ##query the db for the 'PCA' table to retrieve distinct values within the 'POPULATION' column. '.all' executes query + returns results as a list.
    superpopulations = db.session.query(PCA.SUPERPOPULATION).distinct().all() ##NEW EDIT
    return render_template("indextrial.html", populations=populations, superpopulations=superpopulations) ##NEW EDIT - SUPERPOPULATIONS



@app.route('/PCAplot', methods=['POST'])
def plot():
    selected_populations = request.form.getlist('populations')
    selected_superpopulations = request.form.getlist('superpopulations')

    traces = []

    # Filter PCA data based on selected populations
    for population in selected_populations:
        pca_data = PCA.query.filter_by(POPULATION=population).all()
        pc1_values = [entry.PC1 for entry in pca_data]
        pc2_values = [entry.PC2 for entry in pca_data]
        color = f'hsv({len(traces) * (360 // len(selected_populations))}, 100%, 100%)'
        trace = go.Scatter(x=pc1_values, y=pc2_values, mode='markers',
                           marker=dict(color=color, size=10),
                           name=population)
        traces.append(trace)

    # Filter PCA data based on selected superpopulations
    for superpopulation in selected_superpopulations:
        # Extract the first three characters (3-letter code)
        superpopulation_code = superpopulation[:3]
        pca_data = PCA.query.filter(func.substr(PCA.SUPERPOPULATION, 1, 3) == superpopulation_code).all()
        pc1_values = [entry.PC1 for entry in pca_data]
        pc2_values = [entry.PC2 for entry in pca_data]
        color = f'hsv({len(traces) * (360 // len(selected_superpopulations))}, 100%, 100%)'
        trace = go.Scatter(x=pc1_values, y=pc2_values, mode='markers',
                           marker=dict(color=color, size=10),
                           name=superpopulation)
        traces.append(trace)

    layout = go.Layout(title=None, xaxis=dict(title='PC1'), yaxis=dict(title='PC2'))
    fig = go.Figure(data=traces, layout=layout)
    plot_json = fig.to_json()

    return render_template('plot.html', plot_json=plot_json)




if __name__ == '__main__':
    with app.app_context():
        db.create_all()  # Creates all database tables defined by the db models. within this block to create tables only when script is run
    app.run(debug=True)

