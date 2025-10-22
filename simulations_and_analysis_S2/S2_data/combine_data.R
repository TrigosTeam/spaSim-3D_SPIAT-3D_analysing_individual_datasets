setwd("~/R/SPIAT-3D_benchmarking/simulations_and_analysis_S2/S2_data")
ME_metric_df_lists3D <- readRDS("ME_metric_df_lists3D.RDS")
ME_metric_df_lists2D <- readRDS("ME_metric_df_lists2D.RDS")
MN_metric_df_lists3D <- readRDS("MN_metric_df_lists3D.RDS")
MN_metric_df_lists2D <- readRDS("MN_metric_df_lists2D.RDS")
RE_metric_df_lists3D <- readRDS("RE_metric_df_lists3D.RDS")
RE_metric_df_lists2D <- readRDS("RE_metric_df_lists2D.RDS")
RN_metric_df_lists3D <- readRDS("RN_metric_df_lists3D.RDS")
RN_metric_df_lists2D <- readRDS("RN_metric_df_lists2D.RDS")
SE_metric_df_lists3D <- readRDS("SE_metric_df_lists3D.RDS")
SE_metric_df_lists2D <- readRDS("SE_metric_df_lists2D.RDS")
SN_metric_df_lists3D <- readRDS("SN_metric_df_lists3D.RDS")
SN_metric_df_lists2D <- readRDS("SN_metric_df_lists2D.RDS")

metric_df_lists3D <- list(
  mixed_ellipsoid = ME_metric_df_lists3D$mixed_ellipsoid,
  mixed_network = MN_metric_df_lists3D$mixed_network,
  ringed_ellipsoid = RE_metric_df_lists3D$ringed_ellipsoid,
  ringed_network = RN_metric_df_lists3D$ringed_network,
  separated_ellipsoid = SE_metric_df_lists3D$separated_ellipsoid,
  separated_network = SN_metric_df_lists3D$separated_network
)

metric_df_lists2D <- list(
  mixed_ellipsoid = ME_metric_df_lists2D$mixed_ellipsoid,
  mixed_network = MN_metric_df_lists2D$mixed_network,
  ringed_ellipsoid = RE_metric_df_lists2D$ringed_ellipsoid,
  ringed_network = RN_metric_df_lists2D$ringed_network,
  separated_ellipsoid = SE_metric_df_lists2D$separated_ellipsoid,
  separated_network = SN_metric_df_lists2D$separated_network
)

saveRDS(metric_df_lists3D, "metric_df_lists3D.RDS")
saveRDS(metric_df_lists2D, "metric_df_lists2D.RDS")
