FROM python:3.6-stretch
# Install missing deps
USER root

# Added because stretch has been moved to archive.debian.org
RUN echo "deb http://archive.debian.org/debian stretch main"  > /etc/apt/sources.list

# Install unix requirements for the docker
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev     \
        libxslt1-dev    \
        libopenblas-dev \
        liblapack-dev   \
		less			\
		build-essential \
		gfortran        \
		fort77          \
		wget            \
		cmake           \
        libflint-2.5.2  \
        libflint-dev    \
        libgmp-dev      \
		yasm            \
		xvfb            \
		xauth           \
		ffmpeg          \
        firefox-esr

ENV USER NOMAD
ENV HOME /home/$USER

RUN useradd -ms "/bin/bash" "$USER"
USER $USER
WORKDIR $HOME

USER root

# Add extra src if needed
RUN mkdir /src
COPY src/ /src/

# Copy python package requirements
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install --upgrade pipenv
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY ./utils /utils
# Make executable
RUN chmod +x /utils/*.sh

# Install sundials
# RUN /utils/install_chrome.sh  # Removed this on 12/01/2022 due to some broken dependencies for chrome

# Clean up lists
RUN rm -rf /var/lib/apt/lists/*

# Install sundials
RUN /utils/install_sundials_5_1.sh

# NOTE: Sundials version 5.1.0 + scikits.odes==2.6.3 work together
# Export environment variables from installation
ENV SUNDIALS_INCLUDEDIR ="${HOME}/sundials-5.1.0/include"
ENV SUNDIALS_LIBDIR ="${HOME}/sundials-5.1.0/lib"
ENV SUNDIALS_INST="${HOME}/sundials-5.1.0"

# Install python interface to sundials
#modify this
RUN pip install 'scikits.odes== 2.6.3'
ENV LD_LIBRARY_PATH="${HOME}/sundials-5.1.0/lib"
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/lib"


# If files are provided install cplex and gurobi
COPY ./solvers /solvers
RUN /utils/install_cplex.sh
# Install gurobi
RUN /utils/install_gurobi.sh
# Remove installers
RUN rm -rf /solvers
# Activation of necessary licenses
RUN /utils/activate_gurobi.sh

# Make the workspace folder that will link the sources
RUN mkdir /NOMAD

COPY .bashrc $HOME
RUN chown "$USER" "$HOME/.bashrc"
#RUN chown -R "$USER" "$HOME/.cache"

#Finalizing installation

USER $USER
RUN mkdir ./work
WORKDIR ./work

# Load your package in development mode on startup
ENTRYPOINT ["/bin/bash", "-c", "pip install --user -e /NOMAD/docker/pytfa && \
pip install --user -e /NOMAD/docker/skimpy &&\
pip install --user -e /NOMAD/docker/NRAplus &&\
                                $0 $*"]
#ENTRYPOINT ["/bin/bash", "-c", "pip install --user -e /nomad[ORACLE] && \
#                                $0 $*"]

CMD /bin/bash
