from flask import Flask, session, render_template, url_for, redirect, flash
import pandas as pd
from dash import Dash, html, dcc, Output, Input, no_update
import plotly.graph_objects as go
import uuid
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, RDConfig
from io import BytesIO
import base64
from extras.forms import UploadForm
import numpy as np
import os
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer


def smiles_to_base64_img(SMILES):
    mol = Chem.MolFromSmiles(SMILES)
    image = Draw.MolToImage(mol)
    buffer = BytesIO()
    image.save(buffer, format="PNG")
    encoded_image = base64.b64encode(buffer.getvalue()).decode("utf-8")
    return "data:image/png;base64, " + encoded_image


def create_app():
    server = Flask(__name__)
    server.config["SECRET_KEY"] = "1234"
    server.config["MAX_CONTENT_LENGTH"] = 5 * 1024  # 5 KB max-limit.

    @server.route('/', methods=["GET", "POST"])
    def index():
        form = UploadForm()
        if form.validate_on_submit():
            df = pd.read_csv(form.file.data, names=["SMILES"])

            df = pd.concat([df, pd.DataFrame(np.zeros((len(df), 2)),
                                             columns=["clogP", "sascore"])],
                           axis=1)
            for idx, SMILES in enumerate(df.SMILES):
                try:
                    mol = Chem.MolFromSmiles(SMILES)
                    df.loc[idx, ["clogP"]] = Descriptors.MolLogP(mol)
                    df.loc[idx, ["sascore"]] = sascorer.calculateScore(mol)
                except:
                    flash(f"Error. Please check line {idx+1} of your file.")
                    return render_template("main.html", form=form)
            session["data"] = df.to_json()
            return redirect(url_for("next"))
        return render_template("main.html", form=form)

    @server.route("/next", methods=["GET", "POST"])
    def next():
        df = pd.read_json(session["data"])

        app = Dash(server=server, name=str(uuid.uuid4()),
                   url_base_pathname="/" + str(uuid.uuid4()) + "/")
        fig = go.Figure(data=[
            go.Scatter(x=df["clogP"], y=df["sascore"], mode="markers")])
        fig.update_traces(hoverinfo="none", hovertemplate=None)
        fig.update_xaxes(title_text="clogP")
        fig.update_yaxes(title_text="Synthetic accessibility score")

        app.layout = html.Div([
            html.H1("Interactive Visualization"),
            html.H4("This is created by Dash running inside Flask."),
            html.A("Go back to last page", href=url_for("index")),
            dcc.Graph(id="my-graph", figure=fig, clear_on_unhover=True,
                      style={"width": "180vh", "height": "80vh",
                             "display": "block", "margin-left": "auto",
                             "margin-right": "auto", }),
            dcc.Tooltip(id="graph-tooltip")])

        @ app.callback(
            Output("graph-tooltip", "show"),
            Output("graph-tooltip", "bbox"),
            Output("graph-tooltip", "children"),
            Input("my-graph", "hoverData")
        )
        def display_hover(hoverData):
            if hoverData is None:
                return False, no_update, no_update

            pt = hoverData["points"][0]
            bbox = pt["bbox"]
            num = pt["pointNumber"]

            df_row = df.iloc[num]
            SMILE = df_row["SMILES"]

            children = [
                html.Div([
                    html.Img(src=smiles_to_base64_img(
                        SMILE), style={"width": "250px"})],
                    style={"width": "250px", "white-space": "normal"})]
            return True, bbox, children
        return app.index()
    return server


server = create_app()
server.app_context().push()
server.run(host="0.0.0.0", port=80)
