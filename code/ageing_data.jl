using CSV, XLSX, OrderedCollections, Cleaner, DataFrames, Statistics

"""
Extracts fold changes from proteomics and metabolomics data.
Inputs:
    omics_path - Path to the omics data
    oxphos_protein_path - Path to a list of proteins involved in oxidative phosphorylation
    oxphos_metab_path - Path to a list of metabolites involved in oxidative phosphorylation
Outputs:
    df_reaction_FC - A DataFrame containing fold changes for the proteins
    df_metabolite_FC - A DataFrame containing fold changes for the metabolites
"""
function extract_FCs(omics_path, oxphos_protein_path)
    # Extract fold changes for proteins
    (reaction_names,gene_map,filtered_FC,df_reaction_FC) = extract_proteomics(omics_path,oxphos_protein_path)

    # Extract fold changes for metabolites
    df_metabolite_FC = extract_metabolomics(omics_path)

    return df_reaction_FC, df_metabolite_FC
end

"""
Extracts proteomics data
Inputs:
    omics_path - Path to the omics data
    oxphos_protein_path - Path to a list of proteins involved in oxidative phosphorylation

Outputs:
    reaction_names - A list of reactions in oxphos
    gene_map - A dictionary, where the keys are the reactions and the values are lists of protein subunits of the enzymes catalysing the respective reactions
    filtered_FC - Fold changes for each reaction, filtered to the proteins present in the data
    df_reaction_FC - A DataFrame containing the reaction names and fold changes. Where more than one protein is involved, the fold change is taken to be the average fold change among the proteins.
"""
function extract_proteomics(omics_path,oxphos_protein_path)
    table = XLSX.readtable(omics_path,"proteomicsDE")
    df = DataFrame(table...) |> reinfer_schema |> DataFrame

    # Extract gene protein relationships
    filepath = oxphos_protein_path
    (reaction_names,gene_map) = get_gene_map(filepath)
    # Extract fold changes for reactions
    filtered_FC = OrderedDict(
        r => filter_fold_changes(df,genes) for (r,genes) in gene_map
    )
    averaged_FC = [mean(df.FC) for df in values(filtered_FC)]
    df_reaction_FC = DataFrame(
        "name" => reaction_names, 
        "FC" => averaged_FC
    )
    return (reaction_names,gene_map,filtered_FC,df_reaction_FC)
end

"""
Extracts metabolomics data
Inputs:
    omics_path - Path to the omics data
    oxphos_metab_path - Path to a list of metabolites involved in oxidative phosphorylation

Outputs:
    df_metabolite_FC - A DataFrame containing the metabolite names and fold changes.
"""
function extract_metabolomics(omics_path)
    table = XLSX.readtable(omics_path,"metabolomicsDE")
    df = DataFrame(table...) |> reinfer_schema |> DataFrame

    metab_names = df.name
    metab_FC = df.FC
    metab_p = df[!,"adj.P.Val"]
    df_metabolite_FC = DataFrame(
        "name" => metab_names,
        "FC" => metab_FC,
        "adjP" => metab_p
    )

    return (df_metabolite_FC)
end


"""
Returns a dictionary, where the keys are the reactions and the values are lists of protein subunits of the enzymes catalysing the respective reactions
"""
function get_gene_map(filepath)
    sheetname = "Protein list"
    table = XLSX.readtable(filepath, sheetname)
    df_genes = DataFrame(table...) |> reinfer_schema |> DataFrame
    reaction_names = unique(df_genes.Reaction)
    gene_map = OrderedDict(r => gene_associations(r,df_genes) for r in reaction_names)
    return (reaction_names,gene_map)
end

"""
Filters out the proteomics data based on the genes in df_genes
"""
gene_associations(r,df_genes) = filter("Reaction" => n -> n == r, df_genes)[!,"Approved symbol"]

"""
Extracts the first entry of an array. Otherwise returns NaN.
"""
extract(a) = (length(a)==1) ? a[1] : NaN

"""
Returns the fold changes for each reaction, filtered to the proteins present in the data
"""
filter_fold_changes(df,genes) = subset(df, :name => ByRow(n -> n in genes))