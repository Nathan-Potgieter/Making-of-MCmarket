sim_asset_market <- function(corr,
                             T = 252,
                             mv_dist = "t",
                             mv_df = 3,
                             left_cop_weight = 0,
                             left_cop_param = 4,
                             marginal_dist = "norm",
                             marginal_dist_model = NULL,         # may want to change to a list
                             ts_model = list()
) {
    #Simulating innovations
    inno <- sim_inno(corr = corr,
                     mv_dist = mv_dist,
                     df_ellip = mv_df,
                     left_cop_param = left_cop_param,
                     left_cop_weight = left_cop_weight,
                     marginal_dist = marginal_dist,
                     marginal_dist_model = marginal_dist_model,
                     T = T)

    #creating a date vector
    start_date <- Sys.Date()
    all_days <- seq(start_date, start_date %m+%
                        lubridate::days( ceiling( T*(1+(3/7))) ), by = 1)
    weekdays <- all_days[!weekdays(all_days) %in% c('Saturday','Sunday')][1:T]

    #Applying sim_garch to each column in simdat

    simdat <- inno %>% map_dfc(~sim_garch(ts_model, .x))

    #Creating final df
    simdat %>%
        mutate(date = weekdays, .before = `Asset_1`) %>%
        gather(key=Asset, value = Return, -date)
}