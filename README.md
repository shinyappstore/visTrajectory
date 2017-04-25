# visTrajectory

visTrajectory is a Shiny App for visualizing latent trajectories and their
uncertainties in multivariate and high-dimensional datasets.

The app can be downloaded from here or accessed as a web-based tool from shinyapp.io servers:
https://nlhuong.shinyapps.io/visTrajectory/

See the video in the video_demo directory for details on how to use the app.

visTrajectory is build on a RStan-based-package, buds
(https://github.com/nlhuong/buds). buds maps high-dimensional data with an
inherent ordering to latent 1D locations along an underlying gradient or
trajectory. visTrajectory generates plots forrecovered 1D coordinates and their
associated uncertainties as well as provides 2 and 3D usualizations of the
trajectory.
