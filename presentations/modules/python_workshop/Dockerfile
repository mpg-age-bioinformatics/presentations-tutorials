# Copyright (c) Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing.
# Distributed under the terms of the Modified BSD License.

# Debian stretch 9.4-slim image available on 21 June 2018.
FROM hub.age.mpg.de/bioinformatics/software@sha256:5a8722cb45144df91e6e00593aeec8c9f28018bb659b1183ac43790fe3d09b14

LABEL maintainer "bioinformatics@age.mpg.de"

USER root

RUN /bin/bash -c 'sed -i '"'"'24 s/PATH\ /PATH\ $home/'"'"' /modules/modulefiles/general/python/2.7.15 && \
    apt-get update && \
    apt-get install less && \
    echo '"'"'alias ll="ls -l"'"'"' > /home/mpiage/.bash_aliases && \
    echo '"'"'alias l="ls -la"'"'"' >> /home/mpiage/.bash_aliases'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install numpy'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install pandas'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install matplotlib'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install lifelines'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install sklearn'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install scipy'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install matplotlib_venn'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install seaborn'

RUN /bin/bash -c 'source /etc/profile.d/modules.sh && \
module load jupyterhub && \
pip3 install xlrd'

RUN cd /home/_mpiage/ && \
git clone https://github.com/mpg-age-bioinformatics/presentations-tutorials.git && \
cd presentations-tutorials && \
git checkout gh-pages && \
cd ../ && \
mv presentations-tutorials/presentations/modules/python_workshop /home/_mpiage/ && \
rm -rf presentations

ENV HOFOL ".bashrc .bash_logout .Rprofile .software_container python_workshop"

USER $NB_USER

ENTRYPOINT /bin/bash -c '\
for f in $HOFOL ; do \
  if [[ ! -e /home/mpiage/${f} ]] ; \
    then cp -vr /home/_mpiage/${f} /home/mpiage/${f}; \
  fi \
done \
&& source /home/mpiage/.bashrc \
&& /bin/bash'
WORKDIR /home/mpiage



