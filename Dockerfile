# Starting from base miniconda image
FROM continuumio/miniconda3

# Setting vanilla WORKDIR in destination
WORKDIR /app

# This is the complete set of libraries required
RUN conda install -c anaconda scipy
RUN conda install -c anaconda pandas
RUN conda install -c anaconda numpy
RUN conda install -c conda-forge matplotlib
RUN conda install -c anaconda ipywidgets
RUN conda install -c conda-forge panel
RUN conda install scikit-learn
RUN conda install -c conda-forge statsmodels
RUN conda install -c bokeh ipywidgets_bokeh
RUN pip install clustergrammer2
RUN pip install clustergrammer_widget

# Copy the relevant folder into the container
COPY babylon.py .

# Note: the above copies all the contents (not the folder itself) of dep-test into app folder of container

# Run panel serve to start the app
CMD panel serve --address="0.0.0.0" --port=$PORT babylon.py --allow-websocket-origin=als-clustermap.herokuapp.com
