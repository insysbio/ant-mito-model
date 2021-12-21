
##################### Without HetaSimulator ###############
using Parameters, NLopt, LikelihoodProfiler, CSV, DataFrames

pars = (
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

function ant_model(p)

  _pars = merge(pars,p)
  @unpack K_t_mg, K_th_mg, K_d_mg, K_dh_mg, K_t_h, K_d_h, Mg_o, Mg_i, pH_o, pH_i, T_i, T_o, D_i, D_o, fi, del_D, del_T, k1_ANT, k2_ANT, k3_ANT, k4_ANT, K_D_o_ANT, K_T_o_ANT, A, B, C = _pars
  
  H_o = (exp(-pH_o*log(10))*1e6)
  H_i = (exp(-pH_i*log(10))*1e6)

  T_o_free = (T_o/(1+Mg_o/K_t_mg+H_o/K_t_h+Mg_o*H_o/(K_th_mg*K_t_h)))
  T_i_free = (T_i/(1+Mg_i/K_t_mg+H_i/K_t_h+Mg_i*H_i/(K_th_mg*K_t_h)))
  D_o_free = (D_o/(1+Mg_o/K_d_mg+H_o/K_d_h+Mg_o*H_o/(K_dh_mg*K_d_h)))
  D_i_free = (D_i/(1+Mg_i/K_d_mg+H_i/K_d_h+Mg_i*H_i/(K_dh_mg*K_d_h)))

  k1_ANT_fi = (k1_ANT*exp((A*(-4)+B*(-4)+C)*fi))
  k2_ANT_fi = (k2_ANT*exp((A*(-3)+B*(-4)+C)*fi))
  k3_ANT_fi = (k3_ANT*exp((A*(-4)+B*(-3)+C)*fi))
  k4_ANT_fi = (k4_ANT*exp((A*(-3)+B*(-3)+C)*fi))
  K_D_o_ANT_fi = (K_D_o_ANT*exp(3*del_D*fi))
  K_T_o_ANT_fi = (K_T_o_ANT*exp(4*del_T*fi))
  q_fi = (k3_ANT_fi*K_D_o_ANT_fi*exp(fi)/(k2_ANT_fi*K_T_o_ANT_fi))

  ch = (k1_ANT_fi*T_i_free*T_o_free*q_fi/K_T_o_ANT_fi+k2_ANT_fi*T_i_free*D_o_free*q_fi/K_D_o_ANT_fi+k3_ANT_fi*D_i_free*T_o_free/K_T_o_ANT_fi+k4_ANT_fi*D_i_free*D_o_free/K_D_o_ANT_fi)
  zn = ((T_i_free*q_fi+D_i_free)*(1+T_o_free/K_T_o_ANT_fi+D_o_free/K_D_o_ANT_fi))
  
  ch2 = (k2_ANT_fi*T_i_free*D_o_free*q_fi/K_D_o_ANT_fi+k4_ANT_fi*D_i_free*D_o_free/K_D_o_ANT_fi)
  
  x1 = (T_o/(1+Mg_o/K_t_mg+H_o/K_t_h+Mg_o*H_o/(K_th_mg*K_t_h)))   
  x2 = (T_o/(1+Mg_o/K_d_mg+H_o/K_d_h+Mg_o*H_o/(K_dh_mg*K_d_h)))                                     
  
  ch3 = (k1_ANT_fi*T_i_free*x1*q_fi/K_T_o_ANT_fi+k3_ANT_fi*D_i_free*x1/K_T_o_ANT_fi)
  zn3 = ((T_i_free*q_fi+D_i_free)*(1+x1/K_T_o_ANT_fi+x2/K_D_o_ANT_fi))
  
  Vm_T = (((k1_ANT_fi*T_i_free*q_fi+k3_ANT_fi*D_i_free)/(T_i_free*q_fi+D_i_free))*(1/(1+K_T_o_ANT_fi/K_D_o_ANT_fi)))
  ch4 = (k2_ANT_fi*T_i_free*x2*q_fi/K_D_o_ANT_fi+k4_ANT_fi*D_i_free*x2/K_D_o_ANT_fi)
  #ch_obm = (k2_ANT_fi*T_i_free*D_o_free*q_fi/K_D_o_ANT_fi-k3_ANT_fi*D_i_free*T_o_free/K_T_o_ANT_fi)
  #v_obm_ANT = (ch_obm/zn)
  #ff = (fi*26.5)
  #eff = (ch_obm/ch)

  return (
    v1_ANT = ch/zn,
    v2_ANT = ch2/zn,
    v3_ANT = ch3/zn3,
    v4_ANT = ch4/zn3,
    Vm_T = (((k1_ANT_fi*T_i_free*q_fi+k3_ANT_fi*D_i_free)/(T_i_free*q_fi+D_i_free))*(1/(1+K_T_o_ANT_fi/K_D_o_ANT_fi))),
    Vm_D = ((k2_ANT_fi*T_i_free*q_fi+k4_ANT_fi*D_i_free)/(T_i_free*q_fi+D_i_free))*(1/(1+K_D_o_ANT_fi/K_T_o_ANT_fi)),
  )
end

#sim_ant_model(par, fixed, out, fit=()) = ant_model(merge(fixed,fit))

p1 = [288.1, 162.8, 92.91, 48.34, 27.79]
data1 = [28.74, 25.03, 21.11, 16.47, 12.05]
sim1(t,p=()) = ant_model(merge((fi=0, D_o=0, D_i=0, T_i=10000, T_o=t),p))[:v1_ANT]

p2 = [298.6, 71.66, 32.63, 13.01]
data2 = [8.272, 5.063, 3.209, 1.64]
sim2(t,p=()) = ant_model(merge((fi=6.9, D_o=0, D_i=0, T_i=10000, T_o=t),p))[:v1_ANT]

p3 = [303.7, 159, 83.03, 49.1, 29.02]
data3 = [25.24, 22.46, 17.83, 13.26, 9.769]
sim3(d,p=()) = ant_model(merge((fi=0,D_i=10000,T_o=0,T_i=0,D_o=d),p))[:v1_ANT]

p4 = [271.3, 163.5, 93.89, 51.23, 30.23]
data4 = [23.67, 21.75, 18.47, 14.9, 11.77]
sim4(d,p=()) = ant_model(merge((fi=6.9, D_o=d, D_i=10000, T_i=0, T_o=0),p))[:v1_ANT]

p5 = [100.4, 63.86, 41.15, 25.4, 16.2, 10.21]
data5 = [26.95, 23.19, 19.19, 14.54, 10.83, 8.239]
sim5(d,p=()) = ant_model(merge((fi=6.9, D_o=d, D_i=5000, T_i=5000, T_o=0),p))[:v2_ANT]

p6 = [100.5, 63.14, 40.27, 24.84, 15.57, 9.76]
data6 = [24.72, 19.71, 15.42, 11.48, 8.357, 5.709]
sim6(d,p=()) = ant_model(merge((fi=6.9, D_o=d, D_i=5000, T_i=5000, T_o=100),p))[:v2_ANT]

p7 = [101.3, 62.48, 39.88, 24.64, 15.86, 9.581]
data7 = [21.19, 15.89, 11.89, 8.945, 6.297, 4.061]
sim7(d,p=()) = ant_model(merge((fi=6.9, D_o=d, D_i=5000, T_i=5000, T_o=400),p))[:v2_ANT]

p8  =[98.08, 62.69, 43.41, 25.27, 16.06, 10.26]
data8 = [10.95, 8.945, 7.297, 5.355, 4.061, 3.001]
sim8(d,p=()) = ant_model(merge((fi=0, D_o=d, D_i=5000, T_i=5000, T_o=0),p))[:v2_ANT]

p9 = [99.61, 64.95, 40.07, 24.52, 15.33, 9.467]
data9 = [9.71, 8.533, 5.944, 4.296, 3.119, 2.06]
sim9(d,p=()) = ant_model(merge((fi=0, D_o=d, D_i=5000, T_i=5000, T_o=20),p))[:v2_ANT]

p10 = [97.73, 63.2, 37.93, 25.35, 14.42, 8.128]
data10 = [7.18, 5.179, 3.59, 2.472, 1.589, 0.9416]
sim10(d,p=()) = ant_model(merge((fi=0, D_o=d, D_i=5000, T_i=5000, T_o=100),p))[:v2_ANT]

p11 = [300.5, 164.3, 91.74, 48.9, 29.3]
data11 = [34.55, 30.04, 25.67, 19.78, 14.69]
sim11(t,p=()) = ant_model(merge((fi=6.9, D_i=5000, T_i=5000, T_o=t),p))[:v4_ANT]

p12 = [149.5, 44.74, 22.63, 11.17]
data12 = [9.527, 7.127, 5.236, 3.636]
sim12(t,p=()) = ant_model(merge((fi=0, D_i=5000, T_i=5000, T_o=t),p))[:v4_ANT]

p13 = [180, 57.98, 27.78, 11.44]
data13 = [12.22, 9.236, 7.273, 4.073]
sim13(t,p=()) = ant_model(merge((fi=0, D_i=5000, T_i=5000, T_o=t),p))[:v3_ANT]

p14 = [173.2, 69.07, 30.23, 15.31]
data14 = [2.4, 1.818, 1.382, 1.018]
sim14(t,p=()) = ant_model(merge((fi=6.9, D_i=5000, T_i=5000, T_o=t),p))[:v3_ANT]

p15 = [0, 6.9]
data15 = [14.14, 2.53]
sim15(f,p=()) = ant_model(merge((fi=f, D_i=5000, T_i=5000, D_o=0, T_o=400),p))[:Vm_T]

p16 = [0, 6.9]
data16 = [10.57, 39.68]
sim16(f,p=()) = ant_model(merge((fi=f,),p))[:Vm_D]

p_vec = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16]
data_vec = [data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16]
sim_vec = [sim1,sim2,sim3,sim4,sim5,sim6,sim7,sim8,sim9,sim10,sim11,sim12,sim13,sim14,sim15,sim16]

function loss_component(pfit, p, sim, data)
  loss = 0.
    for i in eachindex(p)
      loss += (sim(p[i],pfit) - data[i])^2
    end
  return loss
end

function loss_func(pfit::NamedTuple)
  loss = 0.
  for i in eachindex(p_vec)
    loss += loss_component(pfit, p_vec[i], sim_vec[i], data_vec[i])
  end
  loss
end


function fit_params(
  p;   
  ftol_abs = 0.0,
  ftol_rel = 1e-4, 
  xtol_rel = 0.0,
  xtol_abs = 0.0, 
  fit_alg = :LN_NELDERMEAD,
  maxeval = 10000
)
  obj(pfit::Vector, g::Vector) = loss_func(NamedTuple{Tuple(keys(p))}(pfit))
  
  opt = Opt(fit_alg, length(p))
  opt.min_objective = obj

  opt.ftol_rel = ftol_rel
  opt.ftol_abs = ftol_abs

  opt.xtol_rel = xtol_rel
  opt.xtol_abs = xtol_abs

  opt.maxeval = maxeval

  return NLopt.optimize(opt, collect(values(p)))
end

fitpars = (
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

(minf, minx, ret) = fit_params(fitpars)

ci = [get_interval(
  minx,
  i,
  (p) -> loss_func(NamedTuple{Tuple(keys(fitpars))}(p)),
  :CICO_ONE_PASS,
  theta_bounds = [(-Inf, Inf) for v in minx], #fill((1e-10,1e10), length(p_names[:scn1])),
  scan_bounds=(-100.0, 1e6),
  #local_alg=:LD_MMA,
  scan_tol=1e-5,
  autodiff=false,
  scale =  fill(:direct, length(fitpars)),
  loss_crit = minf + 3.84
  ) for i in 1:length(fitpars)]

plb_1 = [iv.result[1].value for iv in ci]
pub_1 = [iv.result[2].value for iv in ci]

df = DataFrame(param = collect(keys(fitpars)), optim = minx, lb = plb_1, ub = pub_1)