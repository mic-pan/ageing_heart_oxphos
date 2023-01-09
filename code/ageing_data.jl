using CSV, XLSX, OrderedCollections, Cleaner, DataFrames, Statistics

function extract_FCs(omics_path, oxphos_protein_path, oxphos_metab_path)
    # Extract fold changes for proteins
    (reaction_names,gene_map,filtered_FC,df_reaction_FC) = extract_proteomics(omics_path,oxphos_protein_path)

    # Extract fold changes for metabolites
    df_metabolite_FC = extract_metabolomics(omics_path,oxphos_metab_path)

    return df_reaction_FC, df_metabolite_FC
end

function extract_proteomics(omics_path,oxphos_protein_path)
    df = CSV.read(omics_path,DataFrame)
    df_proteins = subset(df, :type => ByRow(n -> n == "proteins"))

    # Extract gene protein relationships
    filepath = oxphos_protein_path
    (reaction_names,gene_map) = get_gene_map(filepath)
    # Extract fold changes for reactions
    filtered_FC = OrderedDict(
        r => filter_fold_changes(df_proteins,genes) for (r,genes) in gene_map
    )
    averaged_FC = [mean(2 .^ df.logFC) for df in values(filtered_FC)]
    df_reaction_FC = DataFrame(
        "name" => reaction_names, 
        "FC" => averaged_FC
    )
    return (reaction_names,gene_map,filtered_FC,df_reaction_FC)
end

function extract_metabolomics(omics_path,oxphos_metab_path)
    df = CSV.read(omics_path,DataFrame)
    df_metabolites = subset(df, :type => ByRow(n -> n == "metabolites"))

    filepath = oxphos_metab_path
    sheetname = "Metabolite list"
    table = XLSX.readtable(filepath, sheetname)
    df_metabolite_names = DataFrame(table...) |> reinfer_schema |> DataFrame
    metab_id = df_metabolite_names[!,"HMDB ID"]
    metab_names = df_metabolite_names[!,"Informal name"]

    metab_FC = [2^extract(df_metabolites[df_metabolites.name .== id, "logFC"]) for id in metab_id]
    metab_p = [extract(df_metabolites[df_metabolites.name .== id, "adjPvalue"]) for id in metab_id]
    df_metabolite_FC = DataFrame(
        "id" => metab_id,
        "name" => metab_names,
        "FC" => metab_FC,
        "adjP" => metab_p
    )

    return (df_metabolite_FC)
end

function get_gene_map(filepath)
    sheetname = "Protein list"
    table = XLSX.readtable(filepath, sheetname)
    df_genes = DataFrame(table...) |> reinfer_schema |> DataFrame
    reaction_names = unique(df_genes.Reaction)
    gene_map = OrderedDict(r => gene_associations(r,df_genes) for r in reaction_names)
    return (reaction_names,gene_map)
end



gene_associations(r,df_genes) = filter("Reaction" => n -> n == r, df_genes)[!,"Approved symbol"]
extract(a) = (length(a)==1) ? a[1] : NaN
filter_fold_changes(df,genes) = subset(df, :name => ByRow(n -> n in genes))