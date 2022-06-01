# Tirocinio

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/DavideFavaro/Tirocinio.jl.svg?branch=master)](https://travis-ci.com/DavideFavaro/Tirocinio.jl)
[![codecov.io](http://codecov.io/github/DavideFavaro/Tirocinio.jl/coverage.svg?branch=master)](http://codecov.io/github/DavideFavaro/Tirocinio.jl?branch=master)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://DavideFavaro.github.io/Tirocinio.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://DavideFavaro.github.io/Tirocinio.jl/dev)
-->

The project has been developed by Davide Favaro as part of a curricular internship at Regione del Veneto, while studying Computer Science at Ca'Foscari university of Venice, Italy, the project was devised and supported by Dr. Gianluca Salogni, company tutor for Regione del Veneto.
The aim is to create a software in pure Julia capable of analyzing the impact of the diffusiion of chemical substances in open areas of Veneto region, Italy, drawing environmental data from measuring stations of the region and satellite images obtained thanks to the Copernicus project.
The software is based on Envifate plugin for QGIS, for the analysis functions, Copernicus API, for data gathering, and NMSIMGIS toolbox, for sound diffusion analysis.
Some of the main scenarios for the analysis are:
- diffusion of chemicals in aquifers;
- diffusion chemicals in lakes;
- diffusion of chemicals in rivers;
- plumes of gasses;
- light pollution;
- noise pollution;
- thermic pollution of rivers;
- sedimentation of pollutants;

Of the aforementined scenarios only the functions for aquifers, lakes, noises, plumes and sediments are working, the remainder are still in development, as are the GUI and the database component.
