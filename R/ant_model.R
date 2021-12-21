library(dMod)
library(dplyr)
library(ggplot2)

x <- Xt()

observables <- c(
  v1_ANT = "(((k1_ANT * exp((A * -4 + B * -4 + C) * fi)) * (T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (T_o / (1 + Mg_o / K_t_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_t_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi))))) / (K_T_o_ANT * exp(4 * del_T * fi)) + ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (D_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi))))) / (K_D_o_ANT * exp(3 * del_D * fi)) + ((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (T_o / (1 + Mg_o / K_t_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_t_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_th_mg * K_t_h)))) / (K_T_o_ANT * exp(4 * del_T * fi)) + ((k4_ANT * exp((A * -3 + B * -3 + C) * fi)) * (D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (D_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h)))) / (K_D_o_ANT * exp(3 * del_D * fi))) / (((T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi)))) + D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (1 + (T_o / (1 + Mg_o / K_t_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_t_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) / (K_T_o_ANT * exp(4 * del_T * fi)) + (D_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) / (K_D_o_ANT * exp(3 * del_D * fi))))",
  v2_ANT = "(((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (D_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi))))) / (K_D_o_ANT * exp(3 * del_D * fi)) + ((k4_ANT * exp((A * -3 + B * -3 + C) * fi)) * (D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (D_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h)))) / (K_D_o_ANT * exp(3 * del_D * fi))) / (((T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi)))) + D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (1 + (T_o / (1 + Mg_o / K_t_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_t_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) / (K_T_o_ANT * exp(4 * del_T * fi)) + (D_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) / (K_D_o_ANT * exp(3 * del_D * fi))))",                                                                                    
  v3_ANT = "(((k1_ANT * exp((A * -4 + B * -4 + C) * fi)) * (T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (T_o / (1 + Mg_o / K_t_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_t_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi))))) / (K_T_o_ANT * exp(4 * del_T * fi)) + ((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (T_o / (1 + Mg_o / K_t_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_t_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_th_mg * K_t_h)))) / (K_T_o_ANT * exp(4 * del_T * fi))) / (((T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi)))) + D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (1 + (T_o / (1 + Mg_o / K_t_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_t_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) / (K_T_o_ANT * exp(4 * del_T * fi)) + (T_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) / (K_D_o_ANT * exp(3 * del_D * fi))))",                                                                                   
  v4_ANT = "(((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (T_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi))))) / (K_D_o_ANT * exp(3 * del_D * fi)) + ((k4_ANT * exp((A * -3 + B * -3 + C) * fi)) * (D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (T_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h)))) / (K_D_o_ANT * exp(3 * del_D * fi))) / (((T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi)))) + D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) * (1 + (T_o / (1 + Mg_o / K_t_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_t_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) / (K_T_o_ANT * exp(4 * del_T * fi)) + (T_o / (1 + Mg_o / K_d_mg + (exp(-pH_o * log(10)) * 1.0e6) / K_d_h + (Mg_o * (exp(-pH_o * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h))) / (K_D_o_ANT * exp(3 * del_D * fi))))",
  Vm_D   = "(((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi)))) + (k4_ANT * exp((A * -3 + B * -3 + C) * fi)) * (D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h)))) / ((T_i / (1 + Mg_i / K_t_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_t_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_th_mg * K_t_h))) * (((k3_ANT * exp((A * -4 + B * -3 + C) * fi)) * (K_D_o_ANT * exp(3 * del_D * fi)) * exp(fi)) / ((k2_ANT * exp((A * -3 + B * -4 + C) * fi)) * (K_T_o_ANT * exp(4 * del_T * fi)))) + D_i / (1 + Mg_i / K_d_mg + (exp(-pH_i * log(10)) * 1.0e6) / K_d_h + (Mg_i * (exp(-pH_i * log(10)) * 1.0e6)) / (K_dh_mg * K_d_h)))) * (1 / (1 + (K_D_o_ANT * exp(3 * del_D * fi)) / (K_T_o_ANT * exp(4 * del_T * fi))))"
)

pars <- c(
    K_t_mg = 101,       # constant of dissociation ATP-Mg {uM}  {1}
    K_th_mg = 10616,    # constant of dissociation ATPH-Mg {uM} {1}
    K_d_mg = 901,       # constant of dissociation ADP-Mg {uM}  {1}
    K_dh_mg = 38168,    # constant of dissociation ADPH-Mg {uM} {1}
    
    K_t_h = 0.323,     # constant of dissociation ATP-H {uM} {1}
    K_d_h = 0.445,     # constant of dissociation ADP-H {uM} {1}

    Mg_o = 0, # 1000,   # external Mg {uM} {2}
    Mg_i = 0, # 1000,   # internal Mg {uM} {2}
    pH_o = 7.2,        # external pH {--} {2}
    pH_i = 7.8,        # internal pH {--} {2}

    T_i = 5000, # 5000, # internal ATP {uM}
    T_o = 2000, # 3000, # external ATP {uM}
    D_i = 5000, # 5000, # internal ADP {uM}
    D_o = 50,          # external ADP {uM}

    fi = 6.9,          # {--} dimentionless potential
    # ! fitted parameters  see article PMID: 16239329
    del_D = 0,         # {--}
    del_T = 0.07,      # {--}
    k1_ANT = 35,       # {1/min}
    k2_ANT = 10.8,     # {1/min}
    k3_ANT = 21,       # {1/min}
    k4_ANT = 29,       # {1/min}
    K_D_o_ANT = 51,    # {uM}
    K_T_o_ANT = 57,    # {uM}
    A = 0.268, # -3.272670e-01,   # -a1 in article {--}
    B = -0.205, # 2.532980e-01,   # -a2 in article {--}
    C = 0.187 # 2.149946e-01,    # a3 in article {--}
)

g <- Y(observables, 
       f = NULL, 
       parameters = names(pars) 
)

times <- c(0.0)

pars1 <- pars
pars1["T_o"] <- 298.6
pars1["fi"] <- 6.9
pars1["D_o"] <- 0
pars1["D_i"] <- 0
pars1["T_i"] <- 10000

out <- (g*x)(times, pars1)
