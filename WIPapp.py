from flask import Flask, render_template, request, url_for, redirect, flash, jsonify, g, send_from_directory
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired
from flask_sqlalchemy import SQLAlchemy
import sqlite3
from sqlalchemy import create_engine, text, func, or_, and_
from sqlalchemy.orm import sessionmaker, aliased
from flask_migrate import Migrate
import pandas as pd
import json
from collections import Counter
from bokeh.io import show, output_file
from bokeh.resources import CDN
from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, TapTool, BoxZoomTool, ResetTool, BoxSelectTool, WheelZoomTool, PanTool, ColumnDataSource, LabelSet, GraphRenderer, EdgesAndLinkedNodes, NodesAndLinkedEdges, LinearColorMapper, BasicTicker, ColorBar
from bokeh.models.graphs import NodesAndLinkedEdges, StaticLayoutProvider
from bokeh.palettes import Category10, Spectral4
from bokeh.plotting import from_networkx, figure
from bokeh.transform import linear_cmap, transform
from bokeh.embed import components, file_html
from bokeh.layouts import row
import networkx as nx
import random
import numpy as np
import math
from pyvis.network import Network
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import os
from copy import deepcopy
from colormaps import rdylgn_r
from tooltips import tip
from flask_bootstrap import Bootstrap


# Initialize Flask application
app = Flask(__name__)
app.config['SECRET_KEY'] = 'CPPKEY'
bootstrap = Bootstrap(app)

# Initialize database
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///cheproadjusted.db'
engine = create_engine('sqlite:///cheproadjusted.db')

# Create session to query database
Session = sessionmaker(bind=engine)
session = Session()

db = SQLAlchemy(app)


# Define table models

class Perturbagen(db.Model):
    __tablename__ = 'Perturbagen'

    name = db.Column(db.String(15), primary_key=True)
    chemspider_id = db.Column(db.Integer)
    action = db.Column(db.String(31))
    synonyms = db.Column(db.String(255))


class Phosphosite(db.Model):
    __tablename__ = 'Phosphosite'

    phosphosite_id = db.Column(db.String(31), primary_key=True)
    uniprot_name = db.Column(db.String(15), nullable=False)
    gene_name = db.Column(db.String(15))
    location = db.Column(db.Integer)
    residue = db.Column(db.String(1))
    quantability = db.Column(db.Float)
    observation = db.Column(db.Integer)
    SITE_7_AA = db.Column(db.Text)
    median_percentile = db.Column(db.Float)
    promiscuity_index = db.Column(db.Float)


class Protein(db.Model):
    __tablename__ = 'Protein'

    uniprot_name = db.Column(db.String(15), primary_key=True)
    kinase_name = db.Column(db.String(255), nullable=False, unique=True)
    expressed_in = db.Column(db.String(31))
    uniprot_id = db.Column(db.String(15), nullable=False)
    gene_name = db.Column(db.String(255))
    gene_synonyms = db.Column(db.String(511))
    description = db.Column(db.Text, nullable=False)
    families = db.Column(db.String(2047))
    length = db.Column(db.Integer, nullable=False)
    sequence = db.Column(db.Text, nullable=False)
    sub_cellular_location = db.Column(db.String(255))
    GO_location = db.Column(db.String(255))
    __table_args__ = (db.UniqueConstraint(
        'kinase_name', name='unq_Protein_kinase_name'),)


class HL60(db.Model):
    __tablename__ = 'HL-60'

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    perturbagen = db.Column(db.Text, db.ForeignKey('Perturbagen.name'))
    phosphosite = db.Column(db.Text)
    cell_line = db.Column(db.Text)
    fold_change = db.Column(db.Float)
    p_value = db.Column(db.Float)
    cv = db.Column(db.Float)
    fc = db.Column(db.Float)
    pval_eb = db.Column(db.Float)
    case = db.Column(db.Text)
    n_runs = db.Column(db.Float)
    n_ctr = db.Column(db.Float)
    n_trt = db.Column(db.Float)
    meansig_ctr = db.Column(db.Float)
    meansig_trt = db.Column(db.Float)
    sid_score = db.Column(db.Float)
    perturbagen_rel = db.relationship(
        'Perturbagen', foreign_keys=[perturbagen])


class KPRelationship(db.Model):
    __tablename__ = 'KP_relationship'

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    database_uniprot_accession = db.Column(db.Text)
    uniprot_name = db.Column(db.Text)
    gene = db.Column(db.Text)
    SITE_7_AA = db.Column(db.Text)
    median_percentile = db.Column(db.Float)
    promiscuity_index = db.Column(db.Float)
    kinase = db.Column(db.String(15))
    kinase_abundance = db.Column(db.Float)
    residue = db.Column(db.Text)
    location = db.Column(db.Text)
    source = db.Column(db.Text)
    cell_line = db.Column(db.Text)
    phosphosite = db.Column(
        db.Text, db.ForeignKey('Phosphosite.phosphosite_id'))
    confidence = db.Column(db.Float)
    phosphosite_id = db.Column(
        db.Text)
    kinase_foreign_key = db.Column(
        db.Text, db.ForeignKey('Protein.kinase_name'))
    kinase_rel = db.relationship('Protein', foreign_keys=[kinase_foreign_key])
    phosphosite_rel = db.relationship(
        'Phosphosite', foreign_keys=[phosphosite])


class MCF7(db.Model):
    __tablename__ = 'MCF-7'

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    perturbagen = db.Column(db.Text, db.ForeignKey('Perturbagen.name'))
    phosphosite = db.Column(db.Text)
    cell_line = db.Column(db.Text)
    fold_change = db.Column(db.Float)
    p_value = db.Column(db.Float)
    cv = db.Column(db.Float)
    fc = db.Column(db.Float)
    pval_eb = db.Column(db.Float)
    case = db.Column(db.Text)
    n_runs = db.Column(db.Float)
    n_ctr = db.Column(db.Float)
    n_trt = db.Column(db.Float)
    meansig_ctr = db.Column(db.Float)
    meansig_trt = db.Column(db.Float)
    sid_score = db.Column(db.Float)
    perturbagen_rel = db.relationship(
        'Perturbagen', foreign_keys=[perturbagen])


class NTERA2CloneD1(db.Model):
    __tablename__ = 'NTERA-2 clone D1'

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    perturbagen = db.Column(db.Text, db.ForeignKey('Perturbagen.name'))
    phosphosite = db.Column(db.Text)
    cell_line = db.Column(db.Text)
    fold_change = db.Column(db.Float)
    p_value = db.Column(db.Float)
    cv = db.Column(db.Float)
    fc = db.Column(db.Float)
    pval_eb = db.Column(db.Float)
    case = db.Column(db.Text)
    n_runs = db.Column(db.Float)
    n_ctr = db.Column(db.Float)
    n_trt = db.Column(db.Float)
    meansig_ctr = db.Column(db.Float)
    meansig_trt = db.Column(db.Float)
    sid_score = db.Column(db.Float)
    perturbagen_rel = db.relationship(
        'Perturbagen', foreign_keys=[perturbagen])


class Observation(db.Model):
    __tablename__ = 'Observation'

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    perturbagen = db.Column(db.String(15), db.ForeignKey('Perturbagen.name'))
    substrate = db.Column(db.String(31))
    cell_line = db.Column(db.String(31))
    fold_change = db.Column(db.Float)
    p_value = db.Column(db.Float)
    cv = db.Column(db.Float)
    perturbagen_rel = db.relationship(
        'Perturbagen', foreign_keys=[perturbagen])


class PKrelationship(db.Model):
    __tablename__ = 'PK_relationship'

    kinase = db.Column(db.String(15), db.ForeignKey(
        'Protein.kinase_name'), primary_key=True)
    perturbagen = db.Column(db.String(15), db.ForeignKey(
        'Perturbagen.name'), primary_key=True)
    source = db.Column(db.String(31), primary_key=True)
    score = db.Column(db.Float)
    kinase_rel = db.relationship('Protein', foreign_keys=[kinase])
    perturbagen_rel = db.relationship(
        'Perturbagen', foreign_keys=[perturbagen])


with app.app_context():
    db.create_all()


# Define search functions

class SearchForm(FlaskForm):
    search_term = StringField('Search', validators=[DataRequired()])
    submit = SubmitField('Submit')


def is_kinase(search_term, session):
    print("hi")
    print(f"**{search_term}**")
    # print(Protein)
    result = session.query(Protein).filter_by(
        kinase_name=search_term).first()
    print(result)
    print(f"Is kinase? {result is not None}")
    return result is not None


def is_protein(search_term, session):
    result = session.query(Protein).filter_by(
        uniprot_name=search_term).first()
    print(f"Is protein? {result is not None}")
    return result is not None


def is_perturbagen(search_term, session):
    result = session.query(Perturbagen).filter(func.lower(
        Perturbagen.name) == func.lower(search_term)).first()
    print(f"Is perturbagen? {result is not None}")
    return result is not None


# HELPER FUNCTIONS #

def linear_normalize(value, min_value, max_value):
    if max_value == min_value:
        return 0  # or any other default value you want to use when the range is zero
    return (value - min_value) / (max_value - min_value)


def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % (int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))


def get_first_match(select_statement, session):
    result = session.execute(select_statement).fetchone()
    return result


def name_to_link(names, route):
    if route == 'protein':
        for i, name in enumerate(names):
            names[i] = '<a href="' + url_for('protein',
                                             protein_name=name) + '">' + name
    elif route == 'kinase':
        for i, name in enumerate(names):
            names[i] = '<a href="' + url_for('kinase',
                                             kinase_name=name) + '">' + name
    return names


# define a simple function to convert p-values to stelar labels
def pstars(p):
    if p is None:
        return ''
    elif p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''


def get_substrates(protein, add_links, session):
    # Construct the SQLAlchemy query using the session.query() method
    # query = session.query(Phosphosite.location, Phosphosite.residue, Phosphosite.phosphosite_id, KPRelationship.location, KPRelationship.residue, KPRelationship.phosphosite_id) \
    # .join(KPRelationship, Phosphosite.phosphosite_id == KPRelationship.phosphosite_id) \
    # .filter(Phosphosite.uniprot_name == protein) \
    # .order_by(Phosphosite.location)

    # *ORIGINAL* Construct the SQLAlchemy query using the session.query() method
    query = session.query(Phosphosite.location, Phosphosite.residue, Phosphosite.phosphosite_id) \
        .filter(Phosphosite.uniprot_name == protein) \
        .order_by(Phosphosite.location)

    # Execute the query
    results = query.all()
    # new_results = []
    # for elem in results:
    # new_results.append((elem[0], elem[1], elem[2]))
    # new_results.append((elem[3], elem[4], elem[5]))
    # Convert query results into a Pandas DF
    results = list(set(results))
    columns = ['location', 'residue', 'phosphosite_id']
    print(results)
    print(columns)
    substrates_df = pd.DataFrame(results, columns=columns)

    detected_in = []
    known_kinases = []
    pdt_kinases = []

    for i, row in substrates_df.iterrows():
        phosphosite_id = row['phosphosite_id']

        # PDTs column
        pdt_query = session.query(KPRelationship.kinase) \
            .filter(KPRelationship.source == 'PDT') \
            .filter(KPRelationship.phosphosite == phosphosite_id)
        pdt_result = pdt_query.all()
        pdt_kinases_list = list(set([r for (r, ) in pdt_result]))
        if add_links:
            pdt_kinases_list = name_to_link(pdt_kinases_list, 'kinase')
        pdt_kinases.append(', '.join(pdt_kinases_list))

        # Reported substrate column
        reported_query = session.query(KPRelationship.kinase) \
            .filter(KPRelationship.source == 'UniProt') \
            .filter(KPRelationship.phosphosite == phosphosite_id)
        reported_result = reported_query.all()
        reported_kinases_list = list(set([r for (r, ) in reported_result]))
        if add_links:
            reported_kinases_list = name_to_link(
                reported_kinases_list, 'kinase')
        known_kinases.append(', '.join(reported_kinases_list))

        # Detected in
        detected_query = session.query(Observation.cell_line) \
            .filter(Observation.substrate == phosphosite_id)
        detected_result = detected_query.all()
        detected_cell_lines = list(set([r for (r, ) in detected_result]))
        cell_line_list = ', '.join(detected_cell_lines)
        cell_line_list = cell_line_list.replace(' clone D1', '')
        detected_in.append(cell_line_list)

    substrates_df['Detected in'] = detected_in
    substrates_df['Reported substrate of'] = known_kinases
    substrates_df['PDT of'] = pdt_kinases

    # Return the DataFrame
    return substrates_df


def kinase_heatmap2(kinase_name, cell_line, DX_THRESHOLD, CV_MAX):
    # retrieve substrate list
    substrates_df = session.query(Phosphosite.phosphosite_id).join(
        KPRelationship, KPRelationship.kinase == kinase_name
    ).filter(
        Phosphosite.phosphosite_id == KPRelationship.phosphosite,
        KPRelationship.source == 'UniProt'
    ).order_by(Phosphosite.phosphosite_id).all()

    substrate_list = [item[0] for item in substrates_df]
    # create fold change matrix to fill with data for specified substrates
    # and as basis of similar p-value and cv matrices
    fc_df = pd.DataFrame()

    for phosphosite_id in substrate_list:
        # Perform the query using the session object and SQLAlchemy models
        query = session.query(Observation.perturbagen, Observation.fold_change, Observation.p_value, Observation.cv).\
            filter(Observation.substrate == phosphosite_id, Observation.cell_line == cell_line).\
            order_by(func.lower(Observation.perturbagen).desc()
                     )  # Using func.lower for case-insensitive sorting

        # Fetch all the query results using the 'all()' method
        results = query.all()

        # Create a DataFrame from the query results
        df = pd.DataFrame(results, columns=[
                          'perturbagen', 'fold_change', 'p_value', 'cv'])
        # print(df)  # Optional: Print the DataFrame to check the results

        if not df.empty:
            if fc_df.empty:
                fc_df['perturbagen'] = df['perturbagen']
                pv_df = fc_df.copy()
                cv_df = fc_df.copy()

            fc_df[phosphosite_id] = df['fold_change']
            pv_df[phosphosite_id] = df['p_value']
            cv_df[phosphosite_id] = df['cv']

    # bail if no observation data available
    if not fc_df.empty:
        fc_df.set_index('perturbagen', inplace=True)
        pv_df.set_index('perturbagen', inplace=True)
        cv_df.set_index('perturbagen', inplace=True)
        # clean up matrix
        H = deepcopy(fc_df)
        # replace sql placeholder with actual values
        H.replace(-888, np.nan, inplace=True)
        H.replace(888, np.inf, inplace=True)
        # replace unreliable (high CV) or cells with NaN
        H.mask(cv_df > CV_MAX, np.nan, inplace=True)
        # create p-value annotations matrix for overlay
        # mask p-values to match fold change dataframe
        pv_df.mask(cv_df > CV_MAX, np.nan, inplace=True)
        pv_df.replace(-888, np.nan, inplace=True)
        pv_annot = pv_df.applymap(pstars)
        # highlight known perturbagens for this kinase
        query = session.query(PKrelationship.perturbagen, PKrelationship.source, PKrelationship.score).\
            filter(PKrelationship.kinase == kinase_name,
                   PKrelationship.score < DX_THRESHOLD)

        # Execute the query and get the result
        result = query.all()

        # 3. Manipulate the data to add the symbol to perturbagen names
        p_list = list(set(row[0] for row in result))
        for i, item in enumerate(H.index.values):
            if item in p_list:
                H.index.values[i] = '\u25BC' + H.index.values[i]

        pv_annot.index = H.index  # make sure p-value index matches
        # create heatmap plot
        # stack for compatibility with bokeh.rect
        fc_stacked_df = H.stack(dropna=False) \
            .rename("fold_change").reset_index()
        pv_stacked_df = pv_annot.stack(dropna=False) \
            .rename("p_value").reset_index()

        mapper = LinearColorMapper(
            palette=rdylgn_r, low=-3, high=3)
        mapper.nan_color = 'lavender'
        # define a figure
        p = figure(
            plot_width=700,
            plot_height=800,
            y_range=list(fc_stacked_df.perturbagen.drop_duplicates()),
            x_range=list(fc_stacked_df.level_1.drop_duplicates()),
            toolbar_location="right",
            tools="save",
            x_axis_location="above")
        # create rectangle for heatmap
        p.rect(
            y="perturbagen",
            x="level_1",
            width=0.95,
            height=0.95,
            source=ColumnDataSource(fc_stacked_df),
            line_color='white',
            fill_color=transform('fold_change', mapper))
        # add legend
        color_bar = ColorBar(
            color_mapper=mapper,
            location=(0, 0),
            ticker=BasicTicker(desired_num_ticks=7),
            major_tick_line_color=None,
            minor_tick_line_color=None)

        p.xaxis.major_label_orientation = math.pi/2
        p.yaxis.major_tick_line_color = None
        p.yaxis.minor_tick_line_color = None
        p.xaxis.axis_line_color = None
        p.yaxis.axis_line_color = None
        p.outline_line_color = None
        p.grid.visible = False
        p.add_layout(color_bar, 'right')
        # add p-stars
        label = LabelSet(y='perturbagen', x='level_1', text='p_value',
                         level='glyph', y_offset=-12,
                         source=ColumnDataSource(pv_stacked_df),
                         text_color='white', text_align='center')
        p.add_layout(label)
    else:
        p = figure(plot_width=700, plot_height=800,
                   toolbar_location=None, tools="")
        p.title.text = ("Insufficient data available for this set of "
                        "perturbagen/substrate combinations")
    # extract component for template and return
    script, div = components(p)
    return {'script': script, 'div': div}


# Define routes

@app.route('/test', methods=['GET', 'POST'])
def test_database():
    try:
        # Perform database operation to test connection
        protein = session.query(Protein).first()

        # If the database operation succeeds, print the retrieved protein
        print(protein)

        return 'Database connection successful'
    except Exception as e:
        print('Database connection error:', str(e))
        return 'Database connection error'


@app.route('/', methods=['GET', 'POST'])
def index():
    form = SearchForm()
    search_term = form.search_term.data
    print(f"Search term: {search_term}")
    if not search_term:
        return render_template('home.html', form=form)
    # return redirect(url_for('kinase', kinase_name=search_term))
    else:
        search_term = search_term.upper()
        if is_kinase(search_term, session):
            return redirect(url_for('kinase', kinase_name=search_term))

        elif is_protein(search_term, session):
            return redirect(url_for('protein', protein_name=search_term))

        elif is_perturbagen(search_term, session):
            return redirect(url_for('perturbagen', perturbagen_name=search_term))
        else:
            flash('Invalid input. Please enter a kinase, protein, or perturbagen name.')
            return render_template('home.html', form=form)


@app.route('/search', methods=['GET', 'POST'])
def search():
    form = SearchForm()
    if form.validate_on_submit():
        search_term = form.search_term.data

        if is_kinase(search_term, session):
            return redirect(url_for('kinase', kinase_name=search_term))

        elif is_protein(search_term, session):
            return redirect(url_for('protein', protein_name=search_term))

        elif is_perturbagen(search_term, session):
            return redirect(url_for('perturbagen', perturbagen_name=search_term))
        else:
            flash('Invalid input. Please enter a kinase, protein, or perturbagen name.')
    return render_template('home.html', form=form)

# want to test table with expandable rows


@app.route('/table')
def table():
    return render_template('table.html')


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/ebdt')
def ebdt():
    return render_template('ebdt.html')


@app.route('/EBDT')
def ebdt_upper():
    return render_template('ebdt.html')


@app.route('/kinase/<kinase_name>', methods=['GET'])
def kinase(kinase_name):
    print(kinase_name)
    kinase = session.query(Protein).filter_by(kinase_name=kinase_name).first()
    print(kinase)
    if kinase is None:
        # Return a 404 error if kinase not found
        return render_template('kinase_error.html'), 404
    else:
        P_MIN = 0.5  # Minimum PDT probability for table
        CV_MAX = 0.5  # Max CV for heatmap
        DX_THRESHOLD = 0.7  # Discoverx significance threshold

        kinase_name = kinase.kinase_name

        query = session.query(
            Phosphosite.phosphosite_id,
            Phosphosite.uniprot_name,
            Phosphosite.residue,
            KPRelationship.source
        ).join(
            KPRelationship,
            Phosphosite.phosphosite_id == KPRelationship.phosphosite_id
        ).filter(
            KPRelationship.kinase == kinase_name
        ).filter(
            or_(KPRelationship.source == 'UniProt',
                KPRelationship.source == 'Cantley')
        )

        results = query.all()

        columns = ['Phosphosite', 'Protein name', 'Residue', 'Source']

        # Step 3: Convert the query result to a DataFrame
        substrates_df = pd.DataFrame(results, columns=columns)

        substrate_count = len(substrates_df)

        # Query the perturbagen data
        perturbagen_data = session.query(Perturbagen.name, Perturbagen.action, PKrelationship.source, PKrelationship.score).\
            join(PKrelationship).\
            filter(PKrelationship.kinase == kinase_name).\
            filter(PKrelationship.score < 0.7).\
            distinct(Perturbagen.name).\
            all()

        merged_perturbagen_data = []
        for perturbagen_name, action, source, score in perturbagen_data:
            vendor = '\u25BC' if source == 'vendor' else ''
            discoverx_assay = int(
                (1 - score) * 100) if source == 'discoverx' else ''
            kuster_et_al = int(
                (1 - score) * 100) if source == 'kuster' else ''

            # Check if the perturbagen already exists in the merged data
            existing_perturbagen = next(
                (p for p in merged_perturbagen_data if p[0]
                 == perturbagen_name), None
            )

            if existing_perturbagen:
                # Update the existing perturbagen entry with additional data
                index = merged_perturbagen_data.index(existing_perturbagen)
                existing_perturbagen[2] = vendor or existing_perturbagen[2]
                existing_perturbagen[3] = discoverx_assay or existing_perturbagen[3]
                existing_perturbagen[4] = kuster_et_al or existing_perturbagen[4]
                merged_perturbagen_data[index] = existing_perturbagen
            else:
                # Add a new entry for the perturbagen
                merged_perturbagen_data.append(
                    [perturbagen_name, action, vendor,
                        discoverx_assay, kuster_et_al]
                )

        perturbagen_count = len(merged_perturbagen_data)

        # Retrieve PDT data from the database
        def retrieve_pdt_data(kinase_query, P_MIN):
            pdt_data = {}
            pdt_count = {}
            graphs = {}
            pdt_svg = {}
            graph_json = {}

            for cell_line in ['MCF-7', 'HL-60', 'NTERA-2 clone D1']:
                query = text(
                    "SELECT Phosphosite.phosphosite_id, Phosphosite.uniprot_name, "
                    "Phosphosite.residue, KP_relationship.confidence "
                    "FROM Phosphosite "
                    "JOIN KP_relationship ON Phosphosite.phosphosite_id = KP_relationship.phosphosite "
                    "JOIN Protein ON Protein.kinase_name = KP_relationship.kinase "
                    "WHERE KP_relationship.kinase = :kinase "
                    "AND KP_relationship.source = 'PDT' "
                    "AND KP_relationship.confidence > :p_min "
                    "AND KP_relationship.cell_line = :cell_line "
                    "ORDER BY Phosphosite.phosphosite_id;"
                )
                cursor = session.execute(
                    query, {"kinase": kinase_query, "p_min": P_MIN, "cell_line": cell_line})

                pdt_df = pd.DataFrame(cursor.fetchall())

                if pdt_df.empty:
                    # Create an empty DataFrame with the desired column names
                    pdt_df = pd.DataFrame(
                        columns=['Phosphosite', 'Protein name', 'Residue', 'Probability'])
                else:
                    # Rename columns if there's data available
                    pdt_df.columns = ['Phosphosite',
                                      'Protein name', 'Residue', 'Probability']

                # Add 'Shared with' column to pdt_df
                pdt_df['Shared with'] = ''

                # check if pdt_df isnt empty before proceeding with loop
                # Find shared kinases for each phosphosite
                connections = []
                for i, row in pdt_df.iterrows():
                    phosphosite = row['Phosphosite']
                    s = text(
                        "SELECT DISTINCT KP_relationship.kinase "
                        "FROM KP_relationship "
                        "WHERE KP_relationship.phosphosite = :phosphosite "
                        "AND KP_relationship.source = 'PDT' "
                        "AND KP_relationship.cell_line = :cell_line "
                        "AND KP_relationship.kinase != :kinase_query"
                    )
                    rp = session.execute(
                        s, {"phosphosite": phosphosite, "cell_line": cell_line, "kinase_query": kinase_query})
                    shared_kinases = [r for (r,) in rp.fetchall()]

                    connections = connections + shared_kinases

                    # Remove duplicates for table
                    shared_kinases = list(set(shared_kinases))
                    shared_kinases = sorted(shared_kinases)
                    shared_with = ', '.join(shared_kinases)
                    pdt_df.at[i, 'Shared with'] = shared_with

                pdt_df[' '] = ''

                pdt_df = pdt_df[[' ', 'Phosphosite', 'Protein name', 'Residue',
                                 'Probability', 'Shared with']]

                pdt_df['Protein name'] = name_to_link(pdt_df['Protein name'],
                                                      'protein')
                pdt_df.rename(columns={' ': tip['vbar']}, inplace=True)

                pdt_data[cell_line] = pdt_df
                pdt_count[cell_line] = len(pdt_df)

                pdt_df.rename(columns={' ': tip['vbar']}, inplace=True)
                pdt_data[cell_line] = pdt_df.to_html(
                    classes='table table-striped table-bordered',
                    table_id=cell_line[0] + 'pdt_table', index=False,
                    escape=False)

                connections = Counter(connections)
                # print(pdt_data)
                # print(connections)
                # print(connections.keys())

                # Graph creation

                # Create a new Pyvis Network object for the current cell line
                net = Network(notebook=False)

                # Add the source kinase (kinase_name) as a node
                # Set the source node color
                net.add_node(kinase_query, color='blue')

                # Define the colormap for the target node colors based on weight
                if connections:
                    min_weight = min(connections.values())
                    max_weight = max(connections.values())
                else:
                    # Set default values for min_weight and max_weight
                    min_weight = 0
                    max_weight = 1
                # RdYlGn is a colormap ranging from red to yellow to green
                cmap = cm.get_cmap('coolwarm')

                # norm = colors.Normalize(vmin=min_weight, vmax=max_weight)

                # Add target kinases as nodes and add edges based on the weight
                for target_kinase, weight in connections.items():
                    # Map the weight to a value between 0 and 1 for the colormap
                    normalized_weight = linear_normalize(
                        weight, min_weight, max_weight)
                    # Get the RGB values from the colormap
                    color = cmap(normalized_weight)[:3]
                    # Convert RGB to hexadecimal format
                    hex_color = rgb_to_hex(color)
                    net.add_node(target_kinase, color=hex_color)

                    net.add_edge(kinase_query, target_kinase, value=weight)

                # net.show_buttons(filter_="physics")

                # Save the Pyvis network graph as an HTML file for visualization
                output_file = os.path.join(
                    "static", f"network_diagram_{cell_line}.html")
                net.save_graph(output_file)

                # Store the Network object itself in the dictionary to access it later
                graphs[cell_line] = net

                # Store the path to the HTML file in the graph data
                # graphs[cell_line] = output_file
            return pdt_data, pdt_count, graphs

        pdt_data, pdt_count, graphs = retrieve_pdt_data(
            kinase_name, P_MIN)

        # print(graphs)

        # get subcell location

        def get_GO_locations(kinase_name):
            # Retrieve GO locations for the specified kinase_name
            go_locations = (
                session.query(Protein.GO_location)
                .filter(Protein.kinase_name == kinase_name)
                .all()
            )

            go_locations = go_locations[0][0]
            go_locations = go_locations.replace(", ", ",")

            # Return the results
            return go_locations

        go_locations = get_GO_locations(kinase_name)

        # Call the function to get GO locations for the specified kinase_name

        # P_MIN = 0.5  # Example value, adjust as needed
        # pdt_data, pdt_count, pdt_svg = retrieve_pdt_data(kinase_name, P_MIN)

        # create the bokeh information for three known substrate heatmaps
        script = {}
        div = {}
        if substrate_count > 0:
            for cell_line in ['MCF-7', 'HL-60', 'NTERA-2 clone D1']:
                heatmap = kinase_heatmap2(kinase_name, cell_line,
                                          DX_THRESHOLD, CV_MAX)
                script[cell_line] = heatmap['script']
                div[cell_line] = heatmap['div']

        return render_template('kinase.html', kinase=kinase, perturbagen_data=merged_perturbagen_data, perturbagen_count=perturbagen_count, graphs=graphs, pdt_data=pdt_data, pdt_count=pdt_count, substrate_count=substrate_count, substrates_df=substrates_df, go_locations=go_locations, script=script, div=div)


@app.route('/phosphosite/vbar/<phosphosite_id>/<cell_line>/<kinase>/<float:DX_THRESHOLD>')
def phosphosite_vbar(phosphosite_id, cell_line, kinase, DX_THRESHOLD):
    # Retrieve observation data
    query = session.query(
        Observation.perturbagen,
        Observation.fold_change,
        Observation.p_value,
        Observation.cv
    ).filter(
        Observation.substrate == phosphosite_id,
        Observation.cell_line == cell_line
    ).order_by(Observation.perturbagen)

    # Execute the query and get the result
    result = query.all()

    # Convert the result to a DataFrame
    columns = ['perturbagen', 'fold_change', 'p_value', 'cv']
    df = pd.DataFrame(result, columns=columns)

    # Replace values that are not SQL-friendly
    df.replace(-888, np.nan, inplace=True)
    df.replace(888, np.inf, inplace=True)

    # Replace p-value column with p-value stars
    df['p_value'] = df['p_value'].apply(pstars)

    # If a kinase is specified, highlight known perturbagens for this kinase
    if kinase:
        # Use aliased to join PKRelationship with the same table (for self-join)
        PKRelationship_alias = aliased(PKrelationship)

        # Perform the query using the session object and the PKRelationship model
        query = session.query(PKrelationship.perturbagen).filter(
            PKrelationship.kinase == kinase,
            PKrelationship.score < DX_THRESHOLD
        )

        # Execute the query and get the result
        result = query.distinct().all()

        # Convert the result to a DataFrame
        dfx = pd.DataFrame(result, columns=['perturbagen'])
        p_list = dfx['perturbagen'].drop_duplicates().tolist()
        for i, item in enumerate(df['perturbagen'].values):
            if item in p_list:
                df['perturbagen'][i] = df['perturbagen'][i] + '\u25BC'

    # Create graph
    # Set vertical position of each star - clip to keep in y-axis range
    df['star_pos'] = 0.2 * np.sign(df['fold_change']) + df['fold_change']
    df['star_pos'] = np.clip(df['star_pos'], -2.8, 2.8)

    source = ColumnDataSource(df)

    p = figure(x_range=df['perturbagen'], plot_height=250,
               y_axis_label='Log fold change', plot_width=800, tools="save", toolbar_location="right", y_range=(-3, 3))
    p.yaxis.axis_label_text_font_style = 'normal'  # override default italics
    p.vbar(x='perturbagen', width=0.6, top='fold_change',
           source=source)
    p.xgrid.grid_line_color = None  # don't show vertical grid lines
    p.xaxis.major_label_orientation = "vertical"
    label = LabelSet(x='perturbagen', y='star_pos', text='p_value',
                     level='glyph', source=source, text_align='center',
                     y_offset=-12)
    p.add_layout(label)
    script, div = components(p)
    print(script)

    return render_template('bokeh_plot.html', script=script, div=div)


@app.route('/phosphosite/minimap/<phosphosite_id>/<float:CV_MAX>')
def phosphosite_minimap(phosphosite_id, CV_MAX):
    data = {}
    cell_lines = ['MCF-7', 'HL-60', 'NTERA-2 clone D1']

    for cell_line in cell_lines:
        # Perform the query using the session and the Observation model
        query = session.query(
            Observation.perturbagen,
            Observation.fold_change,
            Observation.p_value,
            Observation.cv
        ).filter(
            Observation.substrate == phosphosite_id,
            Observation.cell_line == cell_line
        ).order_by(Observation.perturbagen)

        # Execute the query and get the result
        result = query.all()

        # Convert the result to a DataFrame
        columns = ['perturbagen', 'fold_change', 'p_value', 'cv']
        df = pd.DataFrame(result, columns=columns)

        # Replace values that are not SQL-friendly
        df.replace(-888, np.nan, inplace=True)
        df.replace(888, np.inf, inplace=True)

        # Replace p-value column with p-value stars
        df['p_value'] = df['p_value'].apply(pstars)

        data[cell_line] = df

        print(data)

    # create fold change dataframe
    fc_df = pd.DataFrame(data={'MCF-7': data['MCF-7']['fold_change'],
                               'HL-60': data['HL-60']['fold_change'],
                               'NTERA-2 clone D1': data['NTERA-2 clone D1']['fold_change']})
    fc_df.index = data['MCF-7']['perturbagen']

    # create p-value dataframe
    pv_df = pd.DataFrame(data={'MCF-7': data['MCF-7']['p_value'],
                               'HL-60': data['HL-60']['p_value'],
                               'NTERA-2 clone D1': data['NTERA-2 clone D1']['p_value']})
    pv_df.index = data['MCF-7']['perturbagen']

    # create CV dataframe
    cv_df = pd.DataFrame(data={'MCF-7': data['MCF-7']['cv'],
                               'HL-60': data['HL-60']['cv'],
                               'NTERA-2 clone D1': data['NTERA-2 clone D1']['cv']})
    cv_df.index = data['MCF-7']['perturbagen']

    # remove poor quality points and zeros
    # create p-value dataframe
    fc_df.mask(cv_df > CV_MAX, np.nan, inplace=True)
    fc_df.mask(fc_df == 0, np.nan, inplace=True)

    if fc_df.empty:
        p = figure(plot_width=800, plot_height=40, tools="save",
                   toolbar_location="right", x_axis_location="above")
        p.title.text = 'No experimental data for ' + phosphosite_id
    else:
        # stack for compatibility with bokeh.rect
        fc_stacked_df = fc_df.stack(dropna=False). \
            rename("fold_change").reset_index()
        pv_stacked_df = pv_df.stack(
            dropna=False).rename("p_value").reset_index()

        mapper = LinearColorMapper(
            palette=rdylgn_r, low=-3, high=3)
        mapper.nan_color = 'lavender'
        # define a figure
        p = figure(plot_width=800, plot_height=180, tools="save",
                   toolbar_location="right", background_fill_alpha=0,
                   x_axis_location="above",
                   x_range=list(
                       fc_stacked_df.perturbagen.drop_duplicates()),
                   y_range=list(fc_stacked_df.level_1.drop_duplicates()))
        # create rectangle for heatmap
        p.rect(x="perturbagen", y="level_1", width=0.95, height=0.95,
               source=ColumnDataSource(fc_stacked_df), line_color='white',
               fill_color=transform('fold_change', mapper))
        # add legend
        color_bar = ColorBar(color_mapper=mapper, location=(0, 0),
                             ticker=BasicTicker(desired_num_ticks=7),
                             major_tick_line_color=None,
                             minor_tick_line_color=None)

        p.xaxis.major_label_orientation = math.pi/2
        p.yaxis.major_tick_line_color = None
        p.yaxis.minor_tick_line_color = None
        p.xaxis.axis_line_color = None
        p.yaxis.axis_line_color = None
        p.outline_line_color = None
        p.grid.visible = False
        p.add_layout(color_bar, 'right')

        # add p-stars
        label = LabelSet(x='perturbagen', y='level_1', text='p_value',
                         level='glyph', source=ColumnDataSource(pv_stacked_df),
                         angle=math.pi/2, x_offset=12, text_color='white',
                         text_align='center')
        p.add_layout(label)

    script, div = components(p)
    return render_template('bokeh_plot.html', script=script, div=div)


@app.route('/protein/<protein_name>', methods=['GET'])
def protein(protein_name):

    # Try to retrieve information from the database using a session
    protein_info = session.query(Protein).filter_by(
        uniprot_name=protein_name).first()
    print(protein_info)

    if not protein_info:   # Protein not in the database, report this to the user
        protein_info = session.query(Protein.uniprot_name).all()
        return render_template('protein_error.html')
    else:   # We have protein data to show
        # Fetch all substrates for this protein
        substrates_df = get_substrates(
            protein_name, add_links=True, session=session)

        # prepare dataframe ready for rendering
        substrates_df.drop(columns=['phosphosite_id'], inplace=True)
        substrates_df.columns = ['Pos', 'Res', 'Detected in',
                                 'Reported substrate of', 'PDT of']
        substrates_df[' '] = ''
        substrates_df = substrates_df[[' ', 'Pos', 'Res', 'Detected in',
                                       'Reported substrate of', 'PDT of']]
        substrates_df.rename(columns={' ': tip['minimap']}, inplace=True)
        pd.set_option('display.max_colwidth', None)
        substrate_data = substrates_df.to_html(
            classes='table table-hover',
            table_id='substrate_table', index=False, escape=False)
        substrate_count = len(substrates_df)

        def get_GO_locations(protein_name):
            # Retrieve GO locations for the specified kinase_name
            go_locations = (
                session.query(Protein.GO_location)
                .filter(Protein.uniprot_name == protein_name)
                .all()
            )

            go_locations = go_locations[0][0]
            go_locations = go_locations.replace(", ", ",")

            # Return the results
            return go_locations

        go_locations = get_GO_locations(protein_name)

        return render_template('protein.html', protein=protein_info,
                               substrate_count=substrate_count, substrate_data=substrate_data, go_locations=go_locations)


@app.route('/protvista_ptms/<uniprot_name>/<uniprot_id>')
def protvista_ptms(uniprot_name, uniprot_id):

    # Get the sequence of the protein because ProtVista needs this for all sources
    s = session.query(Protein.sequence) \
        .filter(Protein.uniprot_id == uniprot_id)

    # Execute the query and fetch the result
    result = session.execute(s).scalar()
    sequence = get_first_match(s, session)[0]

    print(sequence)

    # Retrieve the substrate dataframe for this protein
    substrates_df = get_substrates(
        uniprot_name, add_links=False, session=session)

    print(substrates_df)

    # drop substrates we haven't detected because there alread in protvista
    substrates_df = substrates_df[substrates_df['Detected in'] != '']

    # assign colour of each substrate depending of associated information
    ptm_dict = {'S': 'phosphoserine', 'T': 'phosphothreonine',
                'Y': 'phosphotyrosine'}
    color_list = []
    description_list = []
    for i, row in substrates_df.iterrows():
        if row['PDT of']:
            color_list.append('#FF0000')   # PDTs are red
            description_list.append('Putative downstream target of '
                                    + row['PDT of'])
        else:
            color_list.append('#00FF00')   # detected substrates are green
            description_list.append('Experimentally observed '
                                    + ptm_dict[row['residue']])
    substrates_df['color'] = color_list
    substrates_df['description'] = description_list

    # extact desired columns and reorganise for conversion to JSON
    subset = substrates_df[['location', 'phosphosite_id', 'color',
                            'description']]
    subset.columns = ['begin', 'phosphosite_id', 'color', 'description']
    subset['end'] = subset['begin']
    subset['type'] = 'MOD_RES'
    subset['category'] = 'PTM'
    # reorder columns and subset furhter
    subset = subset[['type', 'category',
                     'begin', 'end', 'color', 'description']]

    features = subset.to_json(orient='records')

    json_string = ('{"sequence": "' + sequence + '","features": '
                   + features + '}')

    return json_string


@app.route('/perturbagen/<perturbagen_name>', methods=['GET'])
def perturbagen(perturbagen_name):
    DX_THRESHOLD = 0.7
    perturbagen = session.query(Perturbagen).filter(
        Perturbagen.name.ilike(perturbagen_name)).first()

    perturbagen_query = perturbagen.name

    # Query the perturbagen data using an explicit join
    result = session.query(
        PKrelationship.kinase,
        PKrelationship.source,
        PKrelationship.score
    ).join(Perturbagen, PKrelationship.perturbagen == Perturbagen.name).filter(
        Perturbagen.name.ilike(perturbagen_query),
        PKrelationship.score < DX_THRESHOLD
    ).all()

    # Process the query results into the desired format
    kinase_data = {}
    for kinase, source, score in result:
        if kinase not in kinase_data:
            kinase_data[kinase] = {'Kinase': kinase, 'Vendor': '',
                                   'DiscoverX assay': '', 'Kuster <i>et al.</i>': ''}
        if source == 'vendor':
            # assume inhibition so down arrow
            kinase_data[kinase]['Vendor'] = '\u25BC'
        elif source == 'discoverx':
            kinase_data[kinase]['DiscoverX assay'] = '{:d}%'.format(
                int(-(1 - score) * 100))
        elif source == 'kuster':
            kinase_data[kinase]['Kuster <i>et al.</i>'] = '{:d}%'.format(
                int(-(1 - score) * 100))

    # Create a list of dictionaries from the merged data
    kinase_list = list(kinase_data.values())
    kinase_count = len(kinase_list)

    # Format synonyms a bit
    if perturbagen.synonyms:  # only do it there are some synonyms
        synonyms = ', '.join(perturbagen.synonyms.split(','))
    else:
        synonyms = None

    # Pass the DataFrame directly to the template
    return render_template('perturbagen.html', perturbagen=perturbagen, kinase_data=kinase_list, synonyms=synonyms, kinase_count=kinase_count)


# helper functions in helpers.py #
if __name__ == '__main__':
    app.run(debug=True)
