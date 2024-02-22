from flask import Flask, render_template, request, jsonify
import pandas as pd
import os
from flask_sqlalchemy import SQLAlchemy
import plotly.graph_objs as go  # Import Plotly library
import sqlite3
import numpy as np
from sklearn.decomposition import PCA
from sqlalchemy import column, inspect, func



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

class ALLELE_FREQUENCY(db.Model):
    __tablename__ = 'ALLELE_FREQUENCY'
    SNP_ID = db.Column(db.String, db.ForeignKey('VARIANTS.SNP_ID'), primary_key=True)
    ACB = db.Column(db.Numeric)
    ASW = db.Column(db.Numeric)
    BEB = db.Column(db.Numeric)
    CDX = db.Column(db.Numeric)
    CEU = db.Column(db.Numeric)
    CHB = db.Column(db.Numeric)
    CHS = db.Column(db.Numeric)
    CLM = db.Column(db.Numeric)
    ESN = db.Column(db.Numeric)
    FIN = db.Column(db.Numeric)
    GBR = db.Column(db.Numeric)
    GIH = db.Column(db.Numeric)
    GWD = db.Column(db.Numeric)
    IBS = db.Column(db.Numeric)
    ITU = db.Column(db.Numeric)
    JPT = db.Column(db.Numeric)
    KHV = db.Column(db.Numeric)
    LWK = db.Column(db.Numeric)
    MSL = db.Column(db.Numeric)
    MXL = db.Column(db.Numeric)
    PEL = db.Column(db.Numeric)
    PJL = db.Column(db.Numeric)
    PUR = db.Column(db.Numeric)
    SIB = db.Column(db.Numeric)
    STU = db.Column(db.Numeric)
    TSI = db.Column(db.Numeric)
    YRI = db.Column(db.Numeric)

class GENOTYPE_FREQUENCIES(db.Model):
    __tablename__ = 'GENOTYPE_FREQUENCIES'
    SNP_ID = db.Column(db.String, db.ForeignKey('VARIANTS.SNP_ID'), primary_key=True)
    ACB = db.Column(db.Numeric)
    ASW = db.Column(db.Numeric)
    BEB = db.Column(db.Numeric)
    CDX = db.Column(db.Numeric)
    CEU = db.Column(db.Numeric)
    CHB = db.Column(db.Numeric)
    CHS = db.Column(db.Numeric)
    CLM = db.Column(db.Numeric)
    ESN = db.Column(db.Numeric)
    FIN = db.Column(db.Numeric)
    GBR = db.Column(db.Numeric)
    GIH = db.Column(db.Numeric)
    GWD = db.Column(db.Numeric)
    IBS = db.Column(db.Numeric)
    ITU = db.Column(db.Numeric)
    JPT = db.Column(db.Numeric)
    KHV = db.Column(db.Numeric)
    LWK = db.Column(db.Numeric)
    MSL = db.Column(db.Numeric)
    MXL = db.Column(db.Numeric)
    PEL = db.Column(db.Numeric)
    PJL = db.Column(db.Numeric)
    PUR = db.Column(db.Numeric)
    SIB = db.Column(db.Numeric)
    STU = db.Column(db.Numeric)
    TSI = db.Column(db.Numeric)
    YRI = db.Column(db.Numeric)


##Flask routes to handle HTTP requests:

@app.route('/') ##'/' route = index pg
def index():
    populations = db.session.query(PCA.POPULATION).distinct().all() ##query the db for the 'PCA' table to retrieve distinct values within the 'POPULATION' column. '.all' executes query + returns results as a list.
    superpopulations = db.session.query(PCA.SUPERPOPULATION).distinct().all() ##NEW EDIT
    columns = ALLELE_FREQUENCY.__table__.columns.keys()[1:]  # Exclude the first column (SNP_ID)
    return render_template("indextrialSNPsearchClinRelev.html", populations=populations, superpopulations=superpopulations, columns=columns)

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


# @app.route('/snp_search', methods=['POST'])
# def snp_search():
#     # Extract values from the form submission
#     search_option = request.form.get('search_option')
#     snp_selected_populations = request.form.getlist('snp_selected_populations')

#     # Perform database queries to retrieve allele frequencies and genotype frequencies
#     allele_frequencies = {}
#     genotype_frequencies = {}
#     for population in snp_selected_populations:
#         allele_frequency_query = ALLELE_FREQUENCY.query.filter_by(SNP_ID=search_option).first()
#         if allele_frequency_query:
#             allele_frequencies[population] = getattr(allele_frequency_query, population)

#         genotype_frequency_query = GENOTYPE_FREQUENCIES.query.filter_by(SNP_ID=search_option).first()
#         if genotype_frequency_query:
#             genotype_frequencies[population] = getattr(genotype_frequency_query, population)

#     # Query for clinical relevance
#     clinical_relevance_query = VARIANTS.query.filter_by(SNP_ID=search_option).first()
#     clinical_relevance = clinical_relevance_query.CLINICAL_RELEVANCE if clinical_relevance_query else None

#     # Format the results
#     results = {
#         'Allele Frequency': allele_frequencies,
#         'Genotype Frequency': genotype_frequencies,
#         'Clinical Relevance': clinical_relevance
#     }

#     # Return the results in JSON format
#     return jsonify(results)



@app.route('/snp_search', methods=['POST'])
def snp_search():
    # Extract form data
    search_option = request.form.get('search_option')
    search_value = request.form.get('search_value')  # Assuming the input fields are named after the search option
    selected_populations = request.form.getlist('snp-populations')  # Retrieve selected populations

    # Determine the column to search in the VARIANTS table
    if search_option == 'SNP ID':
        search_column = VARIANTS.SNP_ID
    elif search_option == 'Genomic Coordinates':
        search_column = VARIANTS.POS
    elif search_option == 'Gene':
        search_column = VARIANTS.GENE
    else:
        # Handle invalid search option
        return 'Invalid search option', 400

    # Query the database to retrieve clinical relevance
    clinical_relevance = db.session.query(VARIANTS.CLINICAL_RELEVANCE).filter(search_column == search_value).first()

    # if clinical_relevance:
    #     return f'Clinical Relevance: {clinical_relevance[0]}'
    # else:
    #     return 'Clinical relevance not found'
    

    if search_option == 'SNP ID':
        allele_frequency_results = ALLELE_FREQUENCY.query.filter_by(SNP_ID=search_value).filter(ALLELE_FREQUENCY.POPULATION_ID.in_(selected_populations)).all()
        genotype_frequency_results = GENOTYPE_FREQUENCIES.query.filter_by(SNP_ID=search_value).filter(GENOTYPE_FREQUENCIES.POPULATION_ID.in_(selected_populations)).all()
    else:
        # Initialize empty lists if the search option is not 'SNP ID'
        allele_frequency_results = []
        genotype_frequency_results = []

    return render_template('snp_search_results.html', clinical_relevance=clinical_relevance, allele_frequency_results=allele_frequency_results, genotype_frequency_results=genotype_frequency_results) 




if __name__ == '__main__':
    with app.app_context():
        db.create_all()  # Creates all database tables defined by the db models. within this block to create tables only when script is run
    app.run(debug=True)

