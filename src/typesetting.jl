export 
  matrix_to_latex, 
  tensor_to_latex_folded

function matrix_to_latex(M)
    rows = [join(row, " & ") for row in eachrow(M)]
    body = join(rows, " \\")
    return "\begin{bmatrix}" * body * "\end{bmatrix}"
end

function tensor_to_latex_folded(G::Array{T,3}; name="G") where T
    d1, d2, d3 = size(G)
    @assert d1 == d2 == d3 "Tensor must be cubic d×d×d"
    rows_latex = String[]
    for i in 1:d1
        # Collect all entries in the i-th “row of slices”
        row_entries = []
        for k in 1:d3
            for j in 1:d2
                push!(row_entries, string(G[i,j,k]))
            end
            if k < d3
                push!(row_entries, "\vrule")  # vertical separator between slices
            end
        end
        # Join entries of this row with &
        push!(rows_latex, join(row_entries, " & "))
    end
    body = join(rows_latex, " \\")
    col_format = join(["c" for _ in 1:(d2*d3 + (d3-1))], "")
    return "\begin{array}{" * col_format * "}" *
           body *
           "\end{array}"
end
