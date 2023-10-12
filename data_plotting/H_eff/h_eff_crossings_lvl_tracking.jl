using DataFrames # this is like pandas
using CSV 
using Plots 
using LaTeXStrings
using LinearAlgebra
using ProgressBars

# Define useful plotting functions

function nearest_cross_to_v_n(df,cross_n,v_n,ϵ_1)
    df = filter(row ->  row.ϵ_1 == ϵ_1, df)

    df_n  = filter(row ->  row.cross_n == cross_n, df)

    sort!(df_n,:Δnn)
    crosses = [df_n.ϵ_1' df_n.ϵ_2']
    distances = dist(v_n, crosses)

    indx = argmin(distances)

    ϵ_1_cross = df_n.ϵ_1[indx]
    ϵ_2_cross = df_n.ϵ_2[indx]
    return Vector([ϵ_1_cross, ϵ_2_cross])
end

function find_cross_ϵ_2(df,cross_n,ϵ_2)
    df = filter(row ->  row.ϵ_2 == ϵ_2, df)

    df_n  = filter(row ->  row.cross_n == cross_n, df)

    sort!(df_n,:Δnn)
    ϵ_1_cross = df_n.ϵ_1[1]
    ϵ_2_cross = df_n.ϵ_2[1]
    return Vector([ϵ_1_cross, ϵ_2_cross])
end

function find_cross_ϵ_1(df,cross_n,ϵ_1)
    df = filter(row ->  row.ϵ_1 == ϵ_1, df)

    df_n  = filter(row ->  row.cross_n == cross_n, df)

    sort!(df_n,:Δnn)
    ϵ_1_cross = df_n.ϵ_1[1]
    ϵ_2_cross = df_n.ϵ_2[1]
    return Vector([ϵ_1_cross, ϵ_2_cross])
end

function dist(v1, v2)
    mat_of_diffs = v1 .- v2
    norms_squared = sum(mat_of_diffs .* mat_of_diffs,dims=2)
    return sqrt.(norms_squared)
end

df = DataFrame(CSV.File("data/h_eff_anticrossings_lvl_tracking.csv"))

N = unique(df.N)
Δ = unique(df.Δ)
K = unique(df.K)

crossings_n = unique(df.cross_n)

ϵ_1_array = unique(df.ϵ_1)
ϵ_2_array = unique(df.ϵ_2)

# reso_n = 4 # number of resonance lines. Dont know how to automatize this yet


# Define crossings data form
dist_cutoff = 1
first_line_array = []
v_n = [3.1,10]
for ϵ_1 in ProgressBar(ϵ_1_array[end:-10:1])
    df_temp = filter(row ->  row.ϵ_1 == ϵ_1, df) #filter at this level for performance

    v_nn = find_cross_ϵ_1(df_temp,1,ϵ_1) #for the FIRST level
    if dist(v_nn,v_n) < dist_cutoff # check new crossing belongs to the same line
        ϵ_1_cross, ϵ_2_cross = v_nn
        v_n = v_nn
        for n in crossings_n
            # the second condition: otherwise it keeps on going until the crossing from the next line dominates
            if ϵ_2_cross < ϵ_2_array[end] 
                df_temp_prime = filter(row -> row.ϵ_2 == ϵ_2_cross, df_temp)
                df_n = filter(row ->  row.cross_n == n, df_temp_prime)
                ΔE_nm = df_n.Δnn[1]
                first_line_array = vcat(first_line_array, [ΔE_nm ϵ_1_cross ϵ_2_cross n 1])
            end

        end
    end
end

# Define crossings data form
dist_cutoff = 2
second_line_array = []
v_n = [6,10]
for ϵ_1 in ProgressBar(ϵ_1_array[end:-10:1])
    df_temp = filter(row ->  row.ϵ_1 == ϵ_1, df) #filter at this level for performance

    v_nn = find_cross_ϵ_1(df_temp,2,ϵ_1) #for the SECOND level
    if dist(v_nn,v_n) < dist_cutoff # check new crossing belongs to the same line
        ϵ_1_cross, ϵ_2_cross = v_nn
        v_n = v_nn
        for n in crossings_n[2:end]
            # the second condition: otherwise it keeps on going until the crossing from the next line dominates
            if ϵ_2_cross < ϵ_2_array[end] 
                df_temp_prime = filter(row -> row.ϵ_2 == ϵ_2_cross, df_temp)
                df_n = filter(row ->  row.cross_n == n, df_temp_prime)
                ΔE_nm = df_n.Δnn[1]
                second_line_array = vcat(second_line_array, [ΔE_nm ϵ_1_cross ϵ_2_cross n 2])
            end

        end
    end
end

df_line_1 = DataFrame(first_line_array, ["gap","ϵ_1","ϵ_2","cross_n","line_n"])
df_line_2 = DataFrame(second_line_array, ["gap","ϵ_1","ϵ_2","cross_n","line_n"])


df_filtered = filter(row ->  row.cross_n == 1, df_line_1)
df_filtered.gap .= df_filtered.gap/maximum(df_filtered.gap)
scatter!(df_filtered.ϵ_1, df_filtered.ϵ_2 , yerror = df_filtered.gap,label = "n_cross="*string(unique(df_filtered.cross_n)[1]),legend = :topleft,markerstrokecolor = :auto,markersize=2,alpha=.5,xlab=L"\epsilon_1",ylab=L"\epsilon_2")#,xlims=(0,12),ylims=(0,10))


df_filtered = filter(row ->  row.cross_n == 4, df_line_1)
# df_filtered.gap .= df_filtered.gap/maximum(df_filtered.gap)
plot!(df_filtered.ϵ_2, df_filtered.gap,label = L"\delta="*string(unique(df_filtered.cross_n)[1]),legend = :topleft,markerstrokecolor = :auto,markersize=2,alpha=.5,xlab=L"\epsilon_2",ylab=L"G",xlims=(0,10),ylims=(0,15),title="Gap crossing,")





# Define crossings data form
lines_array = []
for _=1:reso_n
    for ϵ_2 in ProgressBar(ϵ_2_array)
        for (nn,mm) in level_pairings
            ϵ_1_cross, ϵ_2_cross = find_cross(df,nn,mm,ϵ_2)
            for (n,m) in level_pairings
                ΔE_nm_array = []
                if n >= nn
                    df_temp = filter(row ->  row.ϵ_1 == ϵ_1_cross && row.ϵ_2 == ϵ_2_cross, df)
                    df_n = filter(row ->  row.n == n, df_temp)
                    df_m = filter(row ->  row.n == m, df_temp)
                    ΔE_nm = abs.( df_n.ΔE_n[1] - df_m.ΔE_n[1] )
                    data_array = vcat(data_array, [ΔE_nm ϵ_1_cross ϵ_2_cross n÷2])
                end
            end
        end
    end
end

# Define crossings data form
data_array = []
for ϵ_2 in ProgressBar(ϵ_2_array)
    for (nn,mm) in level_pairings
        ϵ_1_cross, ϵ_2_cross = find_cross_ϵ_2(df,nn,mm,ϵ_2)
        for (n,m) in level_pairings
            ΔE_nm_array = []
            if n >= nn
                df_temp = filter(row ->  row.ϵ_1 == ϵ_1_cross && row.ϵ_2 == ϵ_2_cross, df)
                df_n = filter(row ->  row.n == n, df_temp)
                df_m = filter(row ->  row.n == m, df_temp)
                ΔE_nm = abs.( df_n.ΔE_n[1] - df_m.ΔE_n[1] )
                data_array = vcat(data_array, [ΔE_nm ϵ_1_cross ϵ_2_cross n÷2])
            end
        end
    end
end


df_data = DataFrame(data_array, ["gap","ϵ_1","ϵ_2","cross_n"])

df_filtered = filter(row ->  row.cross_n == 4, df_data)
# df_filtered.gap .= df_filtered.gap/maximum(df_filtered.gap)
scatter(df_filtered.ϵ_1, df_filtered.ϵ_2 , yerror = df_filtered.gap,label = "n_cross="*string(unique(df_filtered.cross_n)[1]),legend = :topleft,markerstrokecolor = :auto,markersize=2,alpha=.5,xlab=L"\epsilon_1",ylab=L"\epsilon_2")#,xlims=(0,12),ylims=(0,10))