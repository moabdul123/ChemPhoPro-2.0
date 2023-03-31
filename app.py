from flask import Flask, render_template, request, url_for
import pandas as pd
import matplotlib
import requests

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('home.html')


@app.route('/kinase')
def kinase():
    return render_template('kinase.html')


@app.route('/protein')
def protein():
    return render_template('protein.html')


@app.route('/perturbagens')
def perturbagens():
    return render_template('perturbagens.html')
