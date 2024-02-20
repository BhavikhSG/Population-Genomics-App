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
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///C:/Users/xavia/OneDrive/Documents/Masters-PG/Bioinformatics MSc/Software Development Group Project/WebApp/instance/DatabaseTest 2.db' ##configures the URI for connecting to the db; specifies path to the SQLite db file
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

##GOT RID OF THIS TABLE:
# class SAMPLE_POP(db.Model):
#     __tablename__ = 'SAMPLE_POP'
#     SAMPLE_ID = db.Column(db.Integer, primary_key=True, unique=True)
#     POPULATION_CODE = db.Column(db.String)
#     SUPERPOPULATION = db.Column(db.String, db.ForeignKey('SUPERPOP.SUPERPOPULATION_ID'))

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

##need to make ALLELE_FREQUENCY and GENOTYPE_FREQUENCIES table.
# class ALLELE_FREQUENCY(db.Model):
#     __tablename__='ALLELE_FREQUENCY'
#     SNP_ID = db.Column(db.String, primary_key=TRUE, db.ForeignKey('VARIANTS.SNP_ID'))



# DATABASE = os.path.join(app.instance_path, 'New DB.db') ##define variable DATABASE. only need to 
    #include this if going to call on this variable 'DATABASE' later on in the application.
##configuring the database connection for the Flask application by providing the path to the SQLite
##database file within the instance folder.
##this db path is then used by the webapp to connect to the SQLite db and query/modify the data


##Flask routes to handle HTTP requests

@app.route('/') ##'/' route = index pg
def index():
    populations = db.session.query(PCA.POPULATION).distinct().all() ##query the db for the 'PCA' table to retrieve distinct values within the 'POPULATION' column. '.all' executes query + returns results as a list.
    superpopulations = db.session.query(PCA.SUPERPOPULATION).distinct().all() ##NEW EDIT
    return render_template("indextrial.html", populations=populations, superpopulations=superpopulations) ##NEW EDIT - SUPERPOPULATIONS



# @app.route('/PCAplot', methods=['POST']) ##'/PCAplot' route handles POST requests to generate PCA plot
# def plot():
#     selected_populations = request.form.getlist('populations')
#     # selected_superpopulations = [superpopulation[0] for superpopulation in selected_superpopulations]
#     selected_superpopulations = request.form.getlist('superpopulations') ##NEW EDIT
#     # selected_superpopulations = [superpopulation[0] for superpopulation in request.form.getlist('superpopulations')]
#     traces = []

#     for i, population in enumerate(selected_populations):
#         print(f"Selected population: {population}") ##
#         pca_data = PCA.query.filter_by(POPULATION=population).all() ##
#         print(f"PCA data for {population}: {pca_data}") ##
#         pca_data = PCA.query.filter_by(POPULATION=population).all() ##queries 'PCA' table to retrieve records where the 'POPULATION' column matches the specified 'population' value(which user selected - form input checkboxes)
#         pc1_values = [entry.PC1 for entry in pca_data]
#         pc2_values = [entry.PC2 for entry in pca_data]
#         color = f'hsv({i * (360 // len(selected_populations))}, 100%, 100%)' ##generate a unique color for each population
#         trace = go.Scatter(x=pc1_values, y=pc2_values, mode='markers', ##create a scatter plot trace for each population
#                         marker=dict(color=color, size=10),
#                         name=population)
#         traces.append(trace)


#     selected_superpopulations = [superpopulation[0] for superpopulation in selected_superpopulations]
#     # selected_superpopulations = [superpopulation for superpopulation in selected_superpopulations if superpopulation[0] in ['AFR', 'EUR', 'EAS', 'SAS', 'AMR']]
#     for j, superpopulation in enumerate(selected_superpopulations):
#         print(f"Selected superpopulation: {superpopulation}") ##
#         pca_data = PCA.query.filter_by(SUPERPOPULATION=superpopulation).all()
#         # pca_data = PCA.query.filter(PCA.SUPERPOPULATION.in_(selected_superpopulations)).all()
#         ## pca_data = PCA.query.all()
#         ## print("All PCA Data:", pca_data)
#         # pca_data = PCA.query.filter(PCA.SUPERPOPULATION == superpopulation).all()
#         # pca_data = PCA.query.filter(PCA.SUPERPOPULATION.in_(selected_superpopulations)).all() ##
#         print(f"PCA data for {superpopulation}: {pca_data}") ##
#         pc1_values = [entry.PC1 for entry in pca_data]
#         pc2_values = [entry.PC2 for entry in pca_data]
#         color = f'hsv({j * (360 // len(selected_superpopulations))}, 100%, 100%)'
#         trace = go.Scatter(x=pc1_values, y=pc2_values, mode='markers',
#                         marker=dict(color=color, size=10),
#                         name=superpopulation)
#         traces.append(trace)

#     layout = go.Layout(title=None, xaxis=dict(title='PC1'), yaxis=dict(title='PC2'))
#     fig = go.Figure(data=traces, layout=layout)
#     plot_json = fig.to_json()

#     return render_template('plot.html', plot_json=plot_json)

# ------------------------------------------------------------------------------






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





















# ----------------------------------------------------------

# @app.route('/PCAplot', methods=['POST'])
# def plot():
#     selected_populations = request.form.getlist('populations')
#     selected_superpopulations = request.form.getlist('superpopulations')
#     print(selected_populations)
#     print(selected_superpopulations)
#     traces = []

#     for population in selected_populations:
#         pca_data = PCA.query.filter_by(POPULATION=population).all()
#         pc1_values = [entry.PC1 for entry in pca_data]
#         pc2_values = [entry.PC2 for entry in pca_data]
#         trace = go.Scatter(x=pc1_values, y=pc2_values, mode='markers', name=population)
#         traces.append(trace)

#     for superpopulation in selected_superpopulations:
#         pca_data = PCA.query.filter_by(SUPERPOPULATION=superpopulation).all()
#         pc1_values = [entry.PC1 for entry in pca_data]
#         pc2_values = [entry.PC2 for entry in pca_data]
#         trace = go.Scatter(x=pc1_values, y=pc2_values, mode='markers', name=superpopulation)
#         traces.append(trace)

#     layout = go.Layout(title=None, xaxis=dict(title='PC1'), yaxis=dict(title='PC2'))
#     fig = go.Figure(data=traces, layout=layout)
#     plot_json = fig.to_json()

#     return render_template('plot.html', plot_json=plot_json)




if __name__ == '__main__':
    with app.app_context():
        db.create_all()  # Creates all database tables defined by the db models. within this block to create tables only when script is run
    app.run(debug=True)

