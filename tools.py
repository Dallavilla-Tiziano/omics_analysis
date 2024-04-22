# Prepare proteome data for analysis
def prepareProteomeDataset(clinical_data_path, proteome_data_path, data_type, confidence_to_drop=['Low', 'Medium'], abundances_start_at=8, output_index='uniprot', max_na_samples=45):
	''' clinical_data_path: path to clinical data
		proteome_data_path: path to proteome data
		data_type: where the protein data come from (PLASMA, PBMC...)
		confidence_to_drop: values of FDR confidence to drop
		abundances_start_at: number of the column where the protein abundances start (counting from 0!)
		output_index: output counts with 'uniprot' or 'gene_symbol' as index 
	'''

	# Accession columns is going to be set as index. Therefore, we need to subtract 1 from abundances_start_at to account for the column set as index
	data_col = abundances_start_at-1
	# Read clinica data
	metaclinical = pd.read_csv(clinical_data_path)

	# Load proteome raw data 
	proteome = pd.read_csv(proteome_data_path, sep='\t')
	# Isoforms are usually indicated by -2,-3 after te uniprot id
	# Isoforms ids are not recognized by the tool that do uniprot/gene symbol conversion
	# Therefore we create a column with the original uniprot id to be able to associate a gene symbol
	proteome['main_uniprot'] = proteome['Accession'].str.split('-', expand=True)[0]
	proteome.set_index('Accession', inplace=True)
	# Exclude from the analysis proteins with low and medium confidence
	for element in confidence_to_drop:
		proteome = proteome[proteome['Protein FDR Confidence: Combined']!=element]
	# PSMs are needed in DE for data correction therefore we save them in  a table
	psms = proteome['# PSMs']
	psms.columns=['PSM']
	# Drop protein with nan abundance in every sample
	proteome.dropna(subset=proteome.columns[data_col:-1], how="all", inplace=True)
	# Drop isoforms (if we consider isoforms) which we couldn't identify precisely (if '# Protein Unique Peptides' is equal to 0)
	# if '# Protein Unique Peptides' is equal to 0 the estimated abundance is set equal for canonical protein and all isoforms
	# in this case we drop all the isoform and consider only the abundance of the canonical form
	proteome.drop_duplicates(subset=proteome.columns[data_col:-1], inplace=True)
	# protein with '# Unique Peptides'=0 have 0 counts in every sample (this step should may have no effect because of the ste where we  Drop protein with nan abundance in every sample)
	proteome = proteome[proteome['# Unique Peptides']!=0]

	# Rename sample columns with the id of the sample
	columns = []
	if data_type=='PBMC':
	    for column in proteome.columns:
	        if 'Abundance' in column:
	            columns.append(column.split(', ')[-1].lower())
	        else:
	            columns.append(column)
	elif data_type=='PLASMA':
	    for column in proteome.columns:
	        if 'Abundance' in column:
	            columns.append(column.split(', ')[-3].lower())
	        else:
	            columns.append(column)
	proteome.columns = columns

	# Check that all samples are in clinical data (This should be modified in the future because it is based on RNAseq clinical data)
	# Loading RNAseq clincal data is useful to be sure to analyse matching samples but separated analysis may benefit from separated clinical data
	samples_found = []
	samples_NOT_found = []
	for i in range(data_col, len(proteome.columns)-1):
	    sample = proteome.columns[i]
	    if sample in list(metaclinical['id_short']):
	        samples_found.append(sample)
	    else:
	        samples_NOT_found.append(sample)
	proteome.drop(samples_NOT_found, axis=1, inplace=True)
	metaclinical_proteome = pd.DataFrame(columns=metaclinical.columns)
	for sample in samples_found:
	    # if not sample in list(metaclinical_proteome['id_short']):
	    metaclinical_proteome = pd.concat([metaclinical_proteome, metaclinical[metaclinical['id_short']==sample]])

	 # Collapsing replicates
	if output_index=='uniprot':
		proteome_collapsed = pd.DataFrame(index=proteome.index)
	else:
		proteome_collapsed = pd.DataFrame(index=proteome['Gene Symbol'])
	for column in proteome.columns[data_col:-1].unique():
	    temp = proteome.loc[:,column]
	    temp.iloc[:,1] = temp.iloc[:,1].fillna(temp.iloc[:,0])
	    temp.iloc[:,0] = temp.iloc[:,0].fillna(temp.iloc[:,1])
	    proteome_collapsed.loc[:,column] = temp.mean(axis=1)

	# Remove proteins with nan in more than max_na_samples
	proteome_collapsed = proteome_collapsed.loc[proteome_collapsed.isna().sum(axis=1)<max_na_samples,:]

	# Drop replicates from clinical data
	metaclinical_proteome.drop_duplicates('id', inplace=True)

	# Return results
	d_res={}
	d_res['clinical'] = metaclinical_proteome
	d_res['proteome_abundance'] = proteome_collapsed
	d_res['psm'] = psms
	return d_res