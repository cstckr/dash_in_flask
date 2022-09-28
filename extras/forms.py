from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import SubmitField


class UploadForm(FlaskForm):
    file = FileField(validators=[FileRequired(), FileAllowed(["txt"])])
    submit = SubmitField("Upload")
