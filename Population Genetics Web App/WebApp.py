from flask import Flask, render_template, request, jsonify
import pandas as pd
import os
from flask_sqlalchemy import SQLAlchemy
import plotly.graph_objs as go  # Import Plotly library
import sqlite3
import numpy as np
from sklearn.decomposition import PCA
from sqlalchemy import func
import matplotlib.pyplot as plt
from io import BytesIO
from sqlalchemy.orm import aliased
from sqlalchemy import text
from flask import session
from itertools import combinations
import seaborn as sns
import base64

app = Flask(__name__)
app.config['SECRET_KEY'] = '78y23hryufsbu!kjnf7&'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///C:/Users/Harde/OneDrive/Desktop/Close to Submission App/NewTest.db' ##configures the URI for connecting to the db; specifies path to the SQLite db file
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SQLALCHEMY_ECHO'] = True



db = SQLAlchemy(app) ##Initialise SQLAlchemy; creates SQLAlchemy object 'db' that is bound to the Flask app; used to define db models + perform db operations


##Database migration: Database models; Defines several classes that inherit from 'db.Model'=base class provided by SQLAlchemy for defining db models. each class=table in db. each attribute=column in table.
class PCA(db.Model):
    __tablename__ = 'PCA'
    SAMPLE_ID = db.Column(db.Integer,db.ForeignKey('SAMPLE_POP.SAMPLE_ID'), primary_key=True)
    PC1 = db.Column(db.Float)
    PC2 = db.Column(db.Float)

class POP_GROUP(db.Model):
    __tablename__ = 'POP_GROUP'
    POPULATION_ID = db.Column(db.String, primary_key=True, unique=True)
    POPULATION_NAME = db.Column(db.String)

class SUPERPOP(db.Model):
    __tablename__ = 'SUPERPOP'
    SUPERPOPULATION_ID = db.Column(db.String, primary_key=True, unique=True)
    SUPERPOPULATION_NAME = db.Column(db.String)

class ADMIXTURE(db.Model):
    __tablename__ = 'ADMIXTURE'
    SAMPLE_ID = db.Column(db.Integer,db.ForeignKey('SAMPLE_POP.SAMPLE_ID'), primary_key=True)
    V1 = db.Column(db.Float)
    V2 = db.Column(db.Float)
    V3 = db.Column(db.Float)
    V4 = db.Column(db.Float)
    V5 = db.Column(db.Float)

class SAMPLE_POP(db.Model):
    __tablename__ = 'SAMPLE_POP'
    SAMPLE_ID = db.Column(db.Integer, db.ForeignKey('ADMIXTURE.SAMPLE_ID'), primary_key=True)
    POPULATION = db.Column(db.String, db.ForeignKey('POP_GROUP.POPULATION_ID'))
    SUPERPOPULATION = db.Column(db.String, db.ForeignKey('SUPERPOP.SUPERPOPULATION_ID'))

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
    populations = db.session.query(SAMPLE_POP.POPULATION).distinct().all() ##query the db for the 'PCA' table to retrieve distinct values within the 'POPULATION' column. '.all' executes query + returns results as a list.
    superpopulations = db.session.query(SAMPLE_POP.SUPERPOPULATION).distinct().all() ##NEW EDIT
    return render_template("indextrial.html", populations=populations, superpopulations=superpopulations) ##NEW EDIT - SUPERPOPULATIONS

####PCA
@app.route('/PCAplot', methods=['POST'])
def plot():
    # Get selected populations and superpopulations from form data
    selected_populations = request.form.getlist('populations')
    selected_superpopulations = request.form.getlist('superpopulations')

    # List to store plot traces
    traces = []

    # Check if at least one population or superpopulation is selected
    if not selected_populations and not selected_superpopulations:
        return "Please select at least one population or superpopulation to visualize", 400

    # Filter PCA data based on selected populations
    for population in selected_populations:
        # Query PCA data and join with SAMPLE_POP without using alias
        pca_data = (
            PCA.query
            .join(SAMPLE_POP, PCA.SAMPLE_ID == SAMPLE_POP.SAMPLE_ID)
            .filter(SAMPLE_POP.POPULATION == population)
            .all()
        )
        pc1_values = [entry.PC1 for entry in pca_data]
        pc2_values = [entry.PC2 for entry in pca_data]
        color = f'hsv({len(traces) * (360 // len(selected_populations))}, 100%, 100%)'
        trace = go.Scatter(x=pc1_values, y=pc2_values, mode='markers',
                           marker=dict(color=color, size=10),
                           name=population)
        traces.append(trace)

    # Filter PCA data based on selected superpopulations
    for superpopulation in selected_superpopulations:
        superpopulation_code = superpopulation[:3]
        # Query PCA data and join with SAMPLE_POP
        pca_data = (
            PCA.query
            .join(SAMPLE_POP, PCA.SAMPLE_ID == SAMPLE_POP.SAMPLE_ID)
            .filter(func.substr(SAMPLE_POP.SUPERPOPULATION, 1, 3) == superpopulation_code)
            .all()
        )
        pc1_values = [entry.PC1 for entry in pca_data]
        pc2_values = [entry.PC2 for entry in pca_data]
        color = f'hsv({len(traces) * (360 // len(selected_superpopulations))}, 100%, 100%)'
        trace = go.Scatter(x=pc1_values, y=pc2_values, mode='markers',
                           marker=dict(color=color, size=10),
                           name=superpopulation)
        traces.append(trace)

    # Define plot layout
    layout = go.Layout(title=None, xaxis=dict(title='PC1'), yaxis=dict(title='PC2'))

    fig = go.Figure(data=traces, layout=layout)
    # Convert the figure to JSON
    plot_json = fig.to_json()

    # Render the template with the plot JSON
    return render_template('plot.html', plot_json=plot_json)






##### ADMIXTURE
#colours of the 5 ancestors for K=5
cluster_colors = {
    'V1': 'blue',
    'V2': 'green',
    'V3': 'red',
    'V4': 'orange',
    'V5': 'purple'
}
###POPULATION QUERY
def query_populations(selected_populations):
    #if no populations selected return none
    if not selected_populations:
        return None
    #dictionary to contain populations and the values for V1 to V5 
    data = {'Population': [], **{f'V{i}': [] for i in range(1, 6)}}
    #iterate through each of the selected population and query the ADMIXTURE(joined with samplepop) data
    for population in selected_populations:
        admixture_data = db.session.query(ADMIXTURE.V1, ADMIXTURE.V2, ADMIXTURE.V3, ADMIXTURE.V4, ADMIXTURE.V5) \
            .join(SAMPLE_POP, ADMIXTURE.SAMPLE_ID == SAMPLE_POP.SAMPLE_ID) \
            .filter(SAMPLE_POP.POPULATION == population) \
            .all()
        
        for entry in admixture_data:
            data['Population'].append(population)
            # Append the values of V1 to V5
            for i, value in enumerate(entry, 1):
                data[f'V{i}'].append(value)
    #dataframe containing the V1-V5 values for the selected populations
    return pd.DataFrame(data)
###SUPERPOP QUERY
def query_superpopulations(selected_superpopulations):
    #if no superpopulations selected return none
    if not selected_superpopulations:
        return None
    #dicitonary to contain superpopulation and its V1-V5 values
    data = {'Superpopulation': [], **{f'V{i}': [] for i in range(1, 6)}}
    #iterate through each of the selected superpopulation and query the ADMIXTURE(joined with samplepop) data
    for superpopulation in selected_superpopulations:
        admixture_data = db.session.query(ADMIXTURE.V1, ADMIXTURE.V2, ADMIXTURE.V3, ADMIXTURE.V4, ADMIXTURE.V5) \
            .join(SAMPLE_POP, ADMIXTURE.SAMPLE_ID == SAMPLE_POP.SAMPLE_ID) \
            .filter(SAMPLE_POP.SUPERPOPULATION == superpopulation) \
            .all()
        
        for entry in admixture_data:
            data['Superpopulation'].append(superpopulation)
            # Append the values of V1 to V5
            for i, value in enumerate(entry, 1):
                data[f'V{i}'].append(value)
    #dataframe containing the V1-V5 values for the selected populations
    return pd.DataFrame(data)

@app.route('/AdmixturePlot', methods=['POST'])
def AdmixturePlot():
    #gets selected populations/superpopulations from user input
    selected_populations = request.form.getlist('populations')
    selected_superpopulations = request.form.getlist('superpopulations')
    #if none are selected return page saying must select at least one
    if not selected_populations and not selected_superpopulations:
        return "Please select at least one population or superpopulation to visualize ADMIXTURE results", 400
    #gets dataframes from functions above depending on input
    population_df = query_populations(selected_populations)
    superpopulation_df = query_superpopulations(selected_superpopulations)

    #creates matplotlib figure
    fig, ax = plt.subplots(figsize=(12, 8))
    #plots barchart for population data
    if population_df is not None:
        groups = population_df.groupby(population_df.columns[0]) #groups data together by name

        for i, (group_name, group) in enumerate(groups):
            bottom = np.zeros(len(group))
            for col in group.columns[1:]:
                ax.bar(i, group[col], bottom=bottom, label=col if i == 0 else "", color=cluster_colors[col], width=0.8)
                bottom += group[col].values
        #sets labels on axis
        ax.set_xlabel('Population', fontsize=14)
        ax.set_ylabel('Ancestry Fraction', fontsize=14)
        ax.set_title('Admixture Analysis Bar Chart', fontsize=16)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([name for name, _ in groups], rotation=45, ha="right")
        ax.legend(title='Ancestors', bbox_to_anchor=(1.05, 1), loc='upper left')
    #plots barchart for superpop data
    elif superpopulation_df is not None:
        groups = superpopulation_df.groupby(superpopulation_df.columns[0])

        for i, (group_name, group) in enumerate(groups):
            bottom = np.zeros(len(group))
            for col in group.columns[1:]:
                ax.bar(i, group[col], bottom=bottom, label=col if i == 0 else "", color=cluster_colors[col], width=0.8)
                bottom += group[col].values
        #sets labels on axis
        ax.set_xlabel('Superpopulation', fontsize=14)
        ax.set_ylabel('Admixture Analysis Bar Chart', fontsize=14)
        ax.set_title('Ancestry Proportions', fontsize=16)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([name for name, _ in groups], rotation=45, ha="right")
        ax.legend(title='Ancestors', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    plot_data = buffer.getvalue()
    buffer.close()

    return plot_data, 200, {'Content-Type': 'image/png'}


####SNP SEARCH
@app.route('/SNPtable', methods=['POST'])
def SNPtable():
    if request.method == 'POST':
        # Get search type, search value, and selected populations from the form
        search_type = request.form.get('search_type')
        search_value = request.form.get('search_value')
        pos_range_start = request.form.get('pos_range_start')
        pos_range_end = request.form.get('pos_range_end')
        selected_populations = request.form.getlist('populations')

        # Determine which filter to apply based on user input
        if search_type == 'snp_id':
            filter_condition = VARIANTS.SNP_ID == search_value
        elif search_type == 'gene':
            filter_condition = VARIANTS.GENE == search_value
        elif search_type == 'pos_range':
            # Add a filter for POS range if both start and end are provided
            if pos_range_start is not None and pos_range_end is not None:
                filter_condition = VARIANTS.POS.between(int(pos_range_start), int(pos_range_end))
            else:
                return 'Invalid POS range', 400
        else:
            # Handle invalid search
            return 'Invalid search', 400


        result = db.session.query(VARIANTS, AlleleFrequency, GenotypeFrequencies). \
            filter(VARIANTS.SNP_ID == AlleleFrequency.SNP_ID). \
            filter(VARIANTS.SNP_ID == GenotypeFrequencies.SNP_ID). \
            filter(filter_condition). \
            add_columns(VARIANTS.POS, VARIANTS.REF, VARIANTS.ALT, VARIANTS.GENE, VARIANTS.CLINICAL_RELEVANCE, * [getattr(AlleleFrequency, population) for population in selected_populations], * [getattr(GenotypeFrequencies, population) for population in selected_populations]). \
            limit(1000).all()


        # Extract relevant information for printing for SNP SEARCH DATAFRAME
        variant_data = [(variant.SNP_ID, variant.POS, variant.REF, variant.ALT, variant.GENE, variant.CLINICAL_RELEVANCE, *allele_frequency_columns) for variant, *allele_frequency_columns in result]




#######FOR PAIRWISE FST MATRIX PLOTTING
##creating FST data table
        FSTDATA = []

        for snp_id, pos, ref, alt, gene, clinical_relevance, *population_columns in variant_data:
            allele_freqs = {pop: freq for pop, freq in zip(selected_populations, population_columns[7:])}
            FSTDATA.append([snp_id, allele_freqs])

# Creating a DataFrame
        df = pd.DataFrame(FSTDATA, columns=['SNP_ID', 'Allele_Frequencies'])

# Process Allele Frequencies
        data = []
        for index, row in df.iterrows():
            snp_id = row['SNP_ID']
            allele_freqs = row['Allele_Frequencies']

            for pop, freqs in allele_freqs.items():
                alt_af, ref_af = map(float, freqs.split(':'))
                data.append([snp_id, alt_af, ref_af, pop])

# Creating a new DataFrame
        new_columns = ['SNP_ID', 'ALT_AF', 'REF_AF', 'POP']
        new_df = pd.DataFrame(data, columns=new_columns)
        pops = new_df['POP'].unique()
        fst_matrix = np.zeros((len(pops), len(pops)))  # Use NumPy for efficiency

# Pre-calculate means for each POP in a dictionary
        pop_means = new_df.groupby('POP')['ALT_AF'].agg('mean')

# 2. Use NumPy for matrix operations
        n_pops = len(pops)
        fst_matrix = np.zeros((n_pops, n_pops))

        for i in range(n_pops):
            for j in range(i+1, n_pops):
                p1 = pop_means[pops[i]]
                p2 = pop_means[pops[j]]
                fst = 1 - (2 * p1 * p2) / ((p1 + p2) ** 2)
                fst_matrix[i, j] = fst
                fst_matrix[j, i] = fst

# Convert to DataFrame
        fst_matrix_df = pd.DataFrame(fst_matrix, index=pops, columns=pops)
        fst_matrix_txt_path = 'static/fst_matrix.txt'
        fst_matrix_df.to_csv(fst_matrix_txt_path, sep='\t', float_format='%.6f')


# PLOTTING FST MATRIX
        print(fst_matrix_df)

# Convert to float 
        fst_matrix_df = fst_matrix_df.astype(float)

# Create heatmap
        fig, ax = plt.subplots()
        sns.heatmap(fst_matrix_df, annot=True, cmap='coolwarm', xticklabels=pops, yticklabels=pops, ax=ax)
        ax.set_title('Pairwise FST Matrix')

# Save the plot
        temp_file_path = 'static/heatmap.png'
        plt.savefig(temp_file_path)

# Pass the image to  HTML template and dataframe search result for SNP search
        return render_template('SNPtable.html', result=result, selected_populations=selected_populations, search_value=search_value, heatmap_path=temp_file_path)
    
    
    # If the request method is GET, render the template without search results
    return render_template('SNPtable.html', result=None, selected_populations=[], search_value=None)


if __name__ == '__main__':
    with app.app_context():
        db.create_all()  # Creates all database tables defined by the db models. within this block to create tables only when script is run
    app.run(debug=True)

