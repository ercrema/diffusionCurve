FROM rocker/r-ver:4.3.3

# Copy the whole project folder
COPY . /home/diffusionCurve

# Install dependency libraries
RUN  . /etc/environment \
 && chmod -R 777 /home/ \
 && apt-get update \
 && apt-get install -y libicu-dev libglpk-dev libxml2-dev pandoc make libssl-dev libgdal-dev gdal-bin libgeos-dev libproj-dev libsqlite3-dev libudunits2-dev libcurl4-openssl-dev --no-install-recommends \
 && R -q -e 'install.packages(c("truncnorm","here","coda","latex2exp","RColoBrewer","dplyr","nimbleCarbon","nimble","sf","rnaturalearth","rcarbon","emdbook"))' \
 && mkdir /home/output

CMD  cd /home/diffusionCurve \
    &&	Rscript data/jp_clean.R \
    &&	Rscript analysis/japan_abot.R \
    &&	Rscript analysis/post_check_jp_abot.R \
    &&  Rscript data/gb_clean.R \
    &&  Rscript analysis/britain_abot.R \
    &&  Rscript analysis/post_check_gb_abot.R \
    &&  Rscript data/burial_clean.R \
    &&  Rscript analysis/burial_icar.R \
    &&  Rscript sim/simulate_sigmoid.R \
    &&  Rscript sim/simulate_icar.R \
    &&  Rscript sim/fit_sim1a.R \
    &&  Rscript sim/fit_sim1b.R \
    &&  Rscript sim/fit_sim2.R \
    &&  Rscript figures_and_tables/figures_main.R \
    &&  Rscript figures_and_tables/figures_esm.R \
    &&  Rscript figures_and_tables/table_main.R \
    &&  mv /home/diffusionCurve/figures_and_tables/* /home/output/ \
    &&  mv /home/diffusionCurve/results/* /home/output/ \
    &&  mv /home/diffusionCurve/sim/results/* /home/output/
    
# STEP 1: Build the Docker Image. Run the following on your terminal

# docker build -t diffusion .

# STEP 2: Run the docker and save the output.

# docker run -v ~/output:/home/output diffusion

# This will generate a folder called output in your home directory containing figures and tables as well as R images containing the results of the case study  and the tactical simulation.

# Please note that the script will require about 120~150 hours for completion on a desktop machine.








