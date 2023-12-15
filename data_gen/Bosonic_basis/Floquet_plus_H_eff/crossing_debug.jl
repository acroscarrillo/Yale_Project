# H_eff 
df_h_eff =  DataFrame(CSV.File("data/h_eff_crossings_comparison.csv"))

ϵ_1_array_h_eff = unique(df_h_eff.ϵ_1)
ϵ_2_array_h_eff = unique(df_h_eff.ϵ_2)

ϵ_1_test_h_eff =  ϵ_1_array[10]
ϵ_2_test_h_eff =  ϵ_2_array[100]

ΔE_temp_h_eff = filter(row ->  row.ϵ_2 ==ϵ_2_test_h_eff && row.ϵ_1 ==ϵ_1_test_h_eff, df_h_eff)

ΔE_temp_h_eff.ΔE_n[2]
ϵ_1_test_h_eff
m_test = ΔE_temp_h_eff.ΔE_n[2]/ϵ_1_test_h_eff
m_theory = 4*√(ϵ_2_test_h_eff)


# Floquet 
df_floq =  DataFrame(CSV.File("data/floquet_crossings_comparison.csv"))

ϵ_1_array_floq = unique(df_floq.ϵ_1)
ϵ_2_array_floq  = unique(df_floq.ϵ_2)

ϵ_1_test_floq  =  ϵ_1_array_floq[10]
ϵ_2_test_floq  =  ϵ_2_array_floq[100]

ΔE_temp_floq = filter(row ->  row.ϵ_2 ==ϵ_2_test_floq && row.ϵ_1 ==ϵ_1_test_floq, df_floq)

ΔE_temp_floq.ΔE_n[2]
ϵ_1_test_floq
m = ΔE_temp_floq.ΔE_n[3]/ϵ_1_test_floq
m_theory = 4*√(ϵ_2_test_floq)



# kissing sanity check
# H_eff 
df_h_eff =  DataFrame(CSV.File("data/h_eff_crossings_comparison.csv"))
ϵ_1_array_h_eff = unique(df_h_eff.ϵ_1)
ϵ_2_array_h_eff = unique(df_h_eff.ϵ_2)

df_temp_h_eff = filter(row -> row.ϵ_1 == 0 && row.ΔE_n < 300, df_h_eff)

scatter(df_temp_h_eff.ϵ_2,df_temp_h_eff.ΔE_n,ms=1,markerstrokewidth=0)

# Floquet
df_floq =  DataFrame(CSV.File("data/floquet_crossings_comparison.csv"))
ϵ_1_array_floq = unique(df_floq.ϵ_1)
ϵ_2_array_floq = unique(df_floq.ϵ_2)

df_temp_floq = filter(row -> row.ϵ_1 == 0 && row.ΔE_n < 300, df_floq)

scatter!(df_temp_floq.ϵ_2[1:end-20],df_temp_floq.ΔE_n[1:end-20],ms=1,markerstrokewidth=0)



# crossings check
# H_eff 
df_h_eff =  DataFrame(CSV.File("data/h_eff_crossings_comparison.csv"))
ϵ_1_array_h_eff = unique(df_h_eff.ϵ_1)
ϵ_2_array_h_eff = unique(df_h_eff.ϵ_2)

df_temp_h_eff = filter(row -> row.ϵ_2 == ϵ_2_array_h_eff[100] && row.ΔE_n < 300, df_h_eff)

scatter(df_temp_h_eff.ϵ_1,df_temp_h_eff.ΔE_n,ms=1,markerstrokewidth=0)

# Floquet
df_floq =  DataFrame(CSV.File("data/floquet_crossings_comparison.csv"))
ϵ_1_array_floq = unique(df_floq.ϵ_1)
ϵ_2_array_floq = unique(df_floq.ϵ_2)

df_temp_floq = filter(row -> row.ϵ_2 == ϵ_2_array_floq[100] && row.ΔE_n < 300, df_floq)

scatter!(df_temp_floq.ϵ_1[1:end-20],df_temp_floq.ΔE_n[1:end-20],ms=1,markerstrokewidth=0)