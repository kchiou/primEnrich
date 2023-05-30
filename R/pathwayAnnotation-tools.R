#!/usr/bin/env Rscript

.getEnsemblVersion = function(version,species) {
	# suppressPackageStartupMessages(require(biomaRt))
	biomaRt::useEnsembl(
		biomart = 'ENSEMBL_MART_ENSEMBL',
		version = gsub('^Ensembl ','',ensembl.version),
		dataset = paste0(species,'_gene_ensembl')
	)
}

#' Look up valid organism codes
#'
#' This function looks up valid organism values for use with specialized databases.
#' The valid codes are found in the `org_code` column of the returned data frame.
#' 
#' @param database Name of the annotation source.
#' @return A data frame containing an `org_code` column and additional metadata.
#' @export
fetchOrganisms = function(
	database=c('ENSEMBL','KEGG','PANTHER','Reactome','DISEASES','DisGeNET','IID')
) {
	database = database[1]
	
	if (toupper(database) == 'ENSEMBL') {
		# suppressPackageStartupMessages(require(biomaRt))
		
		ens = biomaRt::useEnsembl(biomart='ENSEMBL_MART_ENSEMBL')
		org.metadata = biomaRt::listDatasets(ens)
		org.metadata = data.frame(
			org_code = gsub('_gene_ensembl$','',org.metadata$dataset),
			org.metadata
		)
	} else if (toupper(database) == 'KEGG') {
		# suppressPackageStartupMessages(require(KEGGREST))
		
		genome.list = KEGGREST::keggList('genome')
		genome.list = genome.list[grep('^[a-z]+;',genome.list)]
		
		org.metadata = data.frame(
			org_code = gsub('^([a-z]+); .*','\\1',genome.list),
			t_id = names(genome.list),
			description = gsub('^[a-z]+; ','',genome.list))
	} else if (toupper(database) == 'PANTHER') {
		# suppressPackageStartupMessages(require(rbioapi))
		
		org.metadata = rbioapi::rba_panther_info(what='organisms',verbose=FALSE)
		org.metadata = data.frame(
			org_code = org.metadata$short_name,
			org.metadata)
	} else if (toupper(database) == 'REACTOME') {
		# suppressPackageStartupMessages(require(rbioapi))
		
		org.metadata = rbioapi::rba_reactome_species(only_main=TRUE,verbose=FALSE)
		org.metadata = data.frame(
			org_code = org.metadata$abbreviation,
			org.metadata)
	} else if (toupper(database) %in% c('DISEASES','DISGENET')) {		
		org.metadata = data.frame(
			org_code = 'hsapiens',
			short_name = 'human',
			long_name = 'Homo sapiens')
	} else if (toupper(database) == 'IID') {
		org.list = c(human='Homo sapiens',fly='Drosophila melanogaster',mouse='Mus musculus',rat='Rattus norvegicus',worm='Caenorhabditis elegans',yeast='Saccharomyces cerevisiae',alpaca='Vicugna alpaca',chicken='Gallus gallus',cat='Felis catus',cow='Bos taurus',dog='Canis familiaris',duck='Anas platyrhynchos',guinea_pig='Cavia porcellus',horse='Equus caballus',pig='Sus scrofa',rabbit='Oryctolagus cuniculus',sheep='Ovis aries',turkey='Meleagris gallopavo')
		org.metadata = data.frame(
			org_code = names(org.list),
			short_name = names(org.list),
			long_name = as.character(org.list)
		)
	}
	rownames(org.metadata) = NULL
	org.metadata
}

#' Retrieve current versions of annotations from online sources.
#'
#' This function provides specialized workflows for downloading annotations from
#' established databases and datasets. An internet connection is required.
#' 
#' @param database A character string specifying the name of the annotation to retrieve. If multiple values are given, only the first is used. Valid values are `"GO"`, `"KEGG"`, `"PANTHER"`, `"DISEASES"`, `"DisGeNET"`, and `"IID"`.
#' @param ensembl.species A character string specifying the organism for which genes (ENSEMBL) should be annotated. Valid values can be found using `fetchOrganisms("ENSEMBL")$org_code`.
#' @param kegg.species A character string specifying the organism from which to derive KEGG (Kyoto Encyclopedia of Genes and Genomes) annotations. Valid values can be found using `fetchOrganisms("KEGG")$org_code`. Used only when `database` is `"KEGG"`, ignored otherwise.
#' @param panther.species A character string specifying the organism from which to derive PANTHER (Protein ANalysis THrough Evolutionary Relationships) annotations. Valid values can be found using `fetchOrganisms("PANTHER")$org_code`. Used only when `database` is `"PANTHER"`, ignored otherwise.
#' @param reactome.species A character string specifying the organism from which to derive Reactome annotations. Valid values can be found using `fetchOrganisms("Reactome")$org_code`. Used only when `database` is `"Reactome"`, ignored otherwise.
#' @param iid.species A character string specifying the organism from which to derive IID (Integrated Interactions Database) annotations. Valid values can be found using `fetchOrganisms("IID")$org_code`. Used only when `database` is `"IID"`, ignored otherwise.
#' @param go.namespace A character string specifying the namespace for which to limit GO (Gene Ontology) annotations. Must be one or more of `"BP"`, `"CC"`, or `"MF"`. Used only when `database` is `"GO"`, ignored otherwise.
#' @param diseases.confidence A number specifying the confidence score threshold (DISEASES database) for gene-disease associations to include. If `NULL`, no filter is applied. Used only when `database` is `"DISEASES"`, ignored otherwise.
#' @param disgenet.class A character vector specifying one or more disease class(es) (MeSH codes, DisGeNET database) to include. If `"all"` or `NULL`, all classes are included. Used only when `database` is `"DisGeNET"`, ignored otherwise.
#' @param disgenet.source A character vector specifying one or more gene-disease association source(s) (DisGeNET database) to include. If `"all"` or `NULL`, all sources are included. Used only when `database` is `"DisGeNET"`, ignored otherwise.
#' @param disgenet.score A number specifying the score threshold (DisGeNET database) for gene-disease associations to include. If `NULL`, no filter is applied. Used only when `database` is `"DisGeNET"`, ignored otherwise.
#' @param iid.n.methods A number specifying the number threshold (IID database) of methods supporting a protein-protein interaction. If `NULL`, no filter is applied. Used only when `database` is `"IID"`, ignored otherwise.
#' @param iid.n.publications A number specifying the number threshold (IID database) of publications supporting a protein-protein interaction. If `NULL`, no filter is applied. Used only when `database` is `"IID"`, ignored otherwise.
#' @param iid.evidence.type A character vector specifying one or more type(s) of evidence (IID database) to include. Possible values are `"all"`, `"exp"`, `"pred"`, and `"ortho"`. If `"all"` or `NULL`, all evidence types are included. Used only when `database` is `"IID"`, ignored otherwise.
#' @param iid.direction A character vector specifying the direction(s) of evidence (IID database) to include. Possible values are `"all"`, `"-"`, `">"`, `"<"`, and `"><"`. If `"all"` or `NULL`, all directions are included. Used only when `database` is `"IID"`, ignored otherwise.
#' @param ensembl.version A character string specifying the version of the ENSEMBL database from which to derive annotations.
#' @param n.cores A number specifying the number of cores to use.
#' @return An object of pathwayAnnotation class
#' @export
fetchAnnotations = function(
	database=c('GO','KEGG','PANTHER','Reactome','ReactomePPI','DISEASES','DisGeNET','IID'),
	ensembl.species = 'hsapiens',          # Supported codes in fetchOrganisms('ENSEMBL')
	kegg.species = 'hsa',                  # Supported codes in fetchOrganisms('KEGG')
	panther.species = 'HUMAN',             # Supported codes in fetchOrganisms('PANTHER')
	reactome.species = 'HSA',              # Supported codes in fetchOrganisms('Reactome')
	iid.species = 'human',                 # Supported codes in fetchOrganisms('IID')
	go.namespace = NULL,
	diseases.confidence = NULL,
	disgenet.class = NULL,
	disgenet.source = NULL,
	disgenet.score = NULL,
	iid.n.methods = NULL,
	iid.n.publications = NULL,
	iid.evidence.type = NULL,               # Supported codes: "exp", "pred", "ortho"
	iid.direction = NULL,                   # Supported codes: "-" (direction unknown), ">" (target gene interacting with annotated protein), "<" (annotated protein interacting with target gene), "><" (bidirectional)
	ensembl.version = NULL,
	n.cores = NULL
) {
	# suppressPackageStartupMessages(require(biomaRt))
	
	if (is.null(ensembl.version)) {
		ens = biomaRt::useEnsembl(
			biomart = 'ENSEMBL_MART_ENSEMBL',
			dataset=paste0(ensembl.species,'_gene_ensembl')
		)
		ensembl.version = gsub('Ensembl Genes (.+)','Ensembl \\1',subset(biomaRt::listEnsembl(),biomart == 'genes')$version)
	} else {
		ensembl.version = as.character(ensembl.version)
		ens = biomaRt::useEnsembl(
			biomart = 'ENSEMBL_MART_ENSEMBL',
			version = ensembl.version,
			dataset=paste0(ensembl.species,'_gene_ensembl')
		)
		ensembl.version = subset(biomaRt::listEnsemblArchives(),version == ensembl.version)$name
	}

	database = database[1]
	# if (!is.null(n.cores)) suppressPackageStartupMessages(require(parallel))
	
	if (!toupper(database) %in% c('GO','KEGG','PANTHER','REACTOME','REACTOMEPPI','DISEASES','DISGENET','IID')) {
		stop('Argument database must be one of c("GO","KEGG","PANTHER","Reactome","ReactomePPI","DISEASES","DisGeNET","IID")')
	}
	
	if (toupper(database) == 'GO') {
		database = 'GO'
		# suppressPackageStartupMessages(require(GO.db))
		
		ens.go = biomaRt::getBM(
			attributes=c('ensembl_gene_id','go_id'),
			mart = ens)

		gene2go = if (is.null(n.cores)) {
			lapply(
				sort(unique(ens.go$ensembl_gene_id)),
				function(x) {
					out = sort(ens.go[ens.go$ensembl_gene_id == x,'go_id'])
					out[out != '']
				}
			)
		} else {
			parallel::mclapply(
				sort(unique(ens.go$ensembl_gene_id)),
				function(x) {
					out = sort(ens.go[ens.go$ensembl_gene_id == x,'go_id'])
					out[out != '']
				},
				mc.cores=n.cores
			)
		}
		names(gene2go) = sort(unique(ens.go$ensembl_gene_id))
		
		# Clear empty annotations
		gene2go = gene2go[unlist(lapply(gene2go,length)) > 0]
		
		go.ids = sort(unique(unlist(gene2go)))
		go.ids = go.ids[nchar(go.ids) == 10] # Exclude blank GO terms
		go.info = data.frame(
			pathway_id = go.ids,
			pathway_name = AnnotationDbi::Term(go.ids),
			pathway_namespace = AnnotationDbi::Ontology(go.ids),
			stringsAsFactors=FALSE
		)
		
		if (!is.null(go.namespace) && !all(c('BP','MF','CC') %in% go.namespace)) {
			# If restricting to a particular namespace
			go.info = subset(go.info,pathway_namespace %in% go.namespace)
			gene2go.filtered = if (is.null(n.cores)) {
				lapply(gene2go,function(x) intersect(x,go.info$pathway_id))
			} else {
				parallel::mclapply(gene2go,function(x) intersect(x,go.info$pathway_id),mc.cores=n.cores)
			}
			names(gene2go.filtered) = names(gene2go)
			# Clear empty annotations
			gene2go.filtered  = gene2go.filtered[unlist(lapply(gene2go.filtered,length)) > 0]
			gene2go = gene2go.filtered
			
			database = paste0(database,' (',paste(go.namespace,collapse=', '),')')
		}
		
		rownames(go.info) = NULL
		result = createPathwayAnnotationObject(
			database = database,
			organism = ensembl.species,
			annotation = gene2go,
			definition = go.info,
			ensembl.version = ensembl.version,
			annotation.version = ensembl.version
		)
	} else if (toupper(database) == 'KEGG') {
		database = 'KEGG'
		# suppressPackageStartupMessages(require(KEGGREST))
		if (is.null(kegg.species)) stop('kegg.species must not be NULL')
		
		ens.genes = biomaRt::getBM(
			attributes=c('ensembl_gene_id','entrezgene_id'),
			mart = ens)

		kegg.genes = KEGGREST::keggConv(kegg.species,'ncbi-geneid')

		# KEGG gives NCBI Entrez gene IDs so these have to be converted to Ensembl IDs

		entrez.gene.ids = unique(ens.genes$entrezgene_id)
		entrez.gene.ids = entrez.gene.ids[!is.na(entrez.gene.ids)]

		entrez.kegg = intersect(paste0('ncbi-geneid:',entrez.gene.ids),names(kegg.genes))

		kegg.lookup = kegg.genes[entrez.kegg]

		kegg.pathways = KEGGREST::keggLink('pathway',kegg.species)
		# kegg.ko = KEGGREST::keggLink('ko',kegg.species)
		
		kegg.pathway.info = unlist(strsplit(KEGGREST::keggInfo('pathway'),'\\n'))
		kegg.pathway.version = gsub('^path +','',kegg.pathway.info[grep('^path +Release',kegg.pathway.info)])

		kegg.all = data.frame(kegg_gene_id=names(kegg.pathways),kegg_pathway_id=as.character(kegg.pathways),stringsAsFactors=FALSE)

		genes.join = ens.genes
		genes.join$kegg_gene_id = kegg.lookup[paste0('ncbi-geneid:',ens.genes$entrezgene_id)]

		kegg.all = merge(kegg.all,genes.join,by='kegg_gene_id')

		gene2kegg = if (is.null(n.cores)) {
			lapply(split(kegg.all,kegg.all$ensembl_gene_id),function(x) unique(x$kegg_pathway_id))
		} else {
			parallel::mclapply(split(kegg.all,kegg.all$ensembl_gene_id),function(x) unique(x$kegg_pathway_id),mc.cores=n.cores)
		}

		# Put together KEGG definitions
		kegg.info = data.frame(kegg_pathway_id=sort(unique(kegg.all$kegg_pathway_id)),kegg_pathway_name=NA)
		rownames(kegg.info) = kegg.info$kegg_pathway_id

		# The REST API has a limit so we'll break it down into 10 pathways to query at a time
		# to avoid getting kicked off the server

		while (any(is.na(kegg.info$kegg_pathway_name))) {
			these.pathways = rownames(head(subset(kegg.info,is.na(kegg_pathway_name)),10))
			res = KEGGREST::keggGet(these.pathways)
			keys = paste0('path:',unlist(lapply(res,function(x) x$ENTRY)))
			vals = unlist(lapply(res,function(x) x$NAME))
			kegg.info[keys,'kegg_pathway_name'] = vals
			kegg.info[setdiff(these.pathways,keys),'kegg_pathway_name'] = 'UNKNOWN'
		}
		names(kegg.info) = c('pathway_id','pathway_name')
		
		rownames(kegg.info) = NULL
		result = createPathwayAnnotationObject(
			database = database,
			organism = ensembl.species,
			annotation = gene2kegg,
			definition = kegg.info,
			ensembl.version = ensembl.version,
			annotation.version = kegg.pathway.version
		)
	} else if (toupper(database) == 'PANTHER') {
		database = 'PANTHER'
		# suppressPackageStartupMessages(require(rbioapi))
		if (is.null(panther.species)) stop('panther.species must not be NULL')
		
		ens.genes = biomaRt::getBM(
			attributes=c('ensembl_gene_id'),
			mart = ens)
		
		if (!class(panther.species) %in% c('integer','numeric')) {
			organism.lookup = fetchOrganisms('PANTHER')
			panther.species = as.numeric(subset(organism.lookup,org_code == panther.species)$taxon_id)
		}
		
		parse.panther.results = function(x) {
			this.gene = x$accession
			panther.pathways = unlist(lapply(x$annotation_type_list$annotation_data_type[which(unlist(lapply(x$annotation_type_list$annotation_data_type,function(y) if ('content' %in% names(y)) y$content else '')) == 'ANNOT_TYPE_ID_PANTHER_PATHWAY')],function(z) if ('id' %in% names(z$annotation_list$annotation)) z$annotation_list$annotation$id else unlist(lapply(z$annotation_list$annotation,function(a) a$id))))
			if (!is.null(panther.pathways) & length(panther.pathways)) {
				data.frame(gene = this.gene,panther_pathway_id = panther.pathways)
			} else {
				data.frame(gene = character(0L),panther_pathway_id = character(0L))
			}
		}
		
		panther.pathway.version = subset(rbioapi::rba_panther_info('organisms',verbose=FALSE),taxon_id == panther.species)$version
		
		genes.to.do = ens.genes$ensembl_gene_id
		results.list = list()
		i = 1
		while(length(genes.to.do)) {
			these.genes = genes.to.do[1:min(1000,length(genes.to.do))]
			message(paste0('Fetching set of ',min(1000,length(genes.to.do)),' genes. ',length(setdiff(genes.to.do,these.genes)),' genes to go.'))
			these.results = suppressMessages(rbioapi::rba_panther_mapping(these.genes,panther.species))
			these.results.df = do.call(rbind,if (is.null(n.cores)) {
				lapply(these.results$mapped_genes$gene,parse.panther.results)
			} else {
				parallel::mclapply(these.results$mapped_genes$gene,parse.panther.results,mc.cores=n.cores)
			})
			results.list[[i]] = these.results.df
			genes.to.do = setdiff(genes.to.do,these.genes)
			i = i + 1
		}
		panther.all = do.call(rbind,results.list)
		panther.all$panther_gene_id = gsub('^.+?\\|Ensembl=(.+?)\\|.*','\\1',panther.all$gene)

		gene2pthr = if (is.null(n.cores)) {
			lapply(split(panther.all,panther.all$panther_gene_id),function(x) unique(x$panther_pathway_id))
		} else {
			parallel::mclapply(split(panther.all,panther.all$panther_gene_id),function(x) unique(x$panther_pathway_id),mc.cores=n.cores)
		}
		
		panther.info = rbioapi::rba_panther_info('pathways')
		panther.info = subset(panther.info,pathway_id %in% unique(panther.all$panther_pathway_id))

		rownames(panther.info) = NULL
		result = createPathwayAnnotationObject(
			database = database,
			organism = ensembl.species,
			annotation = gene2pthr,
			definition = panther.info,
			ensembl.version = ensembl.version,
			annotation.version = panther.pathway.version
		)
	} else if (toupper(database) == 'REACTOME') {
		database = 'Reactome'
		# suppressPackageStartupMessages(require(rbioapi))
		if (is.null(reactome.species)) stop('reactome.species must not be NULL')
		
		ens.genes = biomaRt::getBM(
			attributes=c('ensembl_gene_id'),
			mart = ens)

		organism.lookup = fetchOrganisms('Reactome')
		if (!class(reactome.species) %in% c('integer','numeric')) {
			reactome.species = subset(organism.lookup,org_code == reactome.species)$displayName
		} else {
			reactome.species = subset(organism.lookup,dbId == reactome.species)$displayName
		}
		
		temp.file = tempfile()
		
		download.file(paste0('https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt'),temp.file,quiet=TRUE,timeout=600)
		reactome.all = read.delim(
			temp.file,
			header=FALSE,
			col.names=c('id','stId','uri','displayName','evidence','taxon'))
		
		reactome.all = subset(reactome.all,taxon == reactome.species)
		
		reactome.pathway.version = paste('Reactome',rbioapi::rba_reactome_version(verbose=FALSE))
		
		reactome.all$reactome_gene_id = reactome.all$id
		names(reactome.all)[2] = 'pathway_id'
		names(reactome.all)[4] = 'pathway_name'
		
		reactome.all$pathway_id = gsub('-',':',reactome.all$pathway_id)

		gene2react = if (is.null(n.cores)) {
			lapply(split(reactome.all,reactome.all$reactome_gene_id),function(x) unique(x$pathway_id))
		} else {
			parallel::mclapply(split(reactome.all,reactome.all$reactome_gene_id),function(x) unique(x$pathway_id),mc.cores=n.cores)
		}
		
		reactome.info = unique(reactome.all[c('pathway_id','pathway_name','uri')])
		reactome.info = subset(reactome.info,pathway_id %in% unique(reactome.all$pathway_id))

		rownames(reactome.info) = NULL
		result = createPathwayAnnotationObject(
			database = database,
			organism = ensembl.species,
			annotation = gene2react,
			definition = reactome.info,
			ensembl.version = ensembl.version,
			annotation.version = reactome.pathway.version
		)
	} else if (toupper(database) == 'DISEASES') {
		database = 'DISEASES'
		if (ensembl.species != 'hsapiens') {
			warning('DISEASES database is only available for human genome. Use liftAnnotation() to convert. Proceeding with organism "hsapiens".')
			
			# Ensembl version has already been set above
			hsap = .getEnsemblVersion(ensembl.version,'hsapiens')
			ensembl.species = 'hsapiens'
		} else {
			hsap = ens
		}
		temp.file = tempfile()
		download.file('https://download.jensenlab.org/human_disease_textmining_filtered.tsv',temp.file,quiet=TRUE,timeout=600)
		diseases.all = read.delim(
			temp.file,
			col.names=c('protein_id','protein_name','disease_id','disease_name','z_score','confidence','disease_uri'),
			header=FALSE)
		diseases.version = paste0('Downloaded ',`attr<-`(as.POSIXct(Sys.time()),'tzone','UTC'),' UTC')
		
		diseases.info = unique(subset(diseases.all,select=c('disease_id','disease_name')))
		diseases.info = diseases.info[order(as.character(unlist(lapply(strsplit(diseases.info$disease_id,':'),function(x) x[1]))),as.integer(unlist(lapply(strsplit(diseases.info$disease_id,':'),function(x) x[2])))),]
		rownames(diseases.info) = NULL
		names(diseases.info) = c('pathway_id','pathway_name')
		
		hsap.genes = biomaRt::getBM(
			attributes=c('ensembl_gene_id','ensembl_peptide_id','external_gene_name'),
			mart = hsap)
			
		diseases.ensembl = subset(diseases.all,grepl('^ENSP[0-9]{11}',protein_id),select=c('protein_id','disease_id','confidence'))
		diseases.ensembl = merge(diseases.ensembl,hsap.genes,by.x='protein_id',by.y='ensembl_peptide_id',all.x=FALSE,all.y=FALSE)
		
		diseases.names = subset(diseases.all,!grepl('^ENSP[0-9]{11}',protein_id),select=c('protein_id','disease_id','confidence'))
		diseases.names = merge(diseases.names,hsap.genes,by.x='protein_id',by.y='external_gene_name',all.x=FALSE,all.y=FALSE)
		
		diseases.all = unique(rbind(
			diseases.ensembl[,c('ensembl_gene_id','disease_id','confidence')],
			diseases.names[,c('ensembl_gene_id','disease_id','confidence')])
		)
		
		if (!is.null(diseases.confidence)) {
			n.before = nrow(unique(diseases.all[c('disease_id','ensembl_gene_id')]))
			diseases.all = subset(diseases.all,confidence >= diseases.confidence)
			n.after = nrow(unique(diseases.all[c('disease_id','ensembl_gene_id')]))
			message(paste0('Removed ',n.before-n.after,' annotations with confidence < ',paste(diseases.confidence,collapse=', ')))
		}
		
		gene2disease = if (is.null(n.cores)) {
			lapply(split(diseases.all,diseases.all$ensembl_gene_id),function(x) unique(x$disease_id))
		} else {
			parallel::mclapply(split(diseases.all,diseases.all$ensembl_gene_id),function(x) unique(x$disease_id),mc.cores=n.cores)
		}
		
		diseases.info = subset(diseases.info,pathway_id %in% unique(diseases.all$disease_id))
		
		rownames(diseases.info) = NULL
		result = createPathwayAnnotationObject(
			database = database,
			organism = ensembl.species,
			annotation = gene2disease,
			definition = diseases.info,
			ensembl.version = ensembl.version,
			annotation.version = diseases.version
		)
	} else if (toupper(database) == toupper('DisGeNET')) {
		database = 'DisGeNET'
		if (ensembl.species != 'hsapiens') {
			warning('DisGeNET database is only available for human genome. Use liftAnnotation() to convert. Proceeding with organism "hsapiens".')
			
			# Ensembl version has already been set above
			hsap = .getEnsemblVersion(ensembl.version,'hsapiens')
			ensembl.species = 'hsapiens'
		} else {
			hsap = ens
		}
		temp.file = tempfile()
		download.file('https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz',temp.file,quiet=TRUE,timeout=600)
		disgenet.all = read.delim(
			gzfile(temp.file),
			header=TRUE)
		disgenet.version = paste0('Downloaded ',`attr<-`(as.POSIXct(Sys.time()),'tzone','UTC'),' UTC')
		
		disgenet.info = unique(disgenet.all[c('diseaseId','diseaseName','diseaseClass','diseaseSemanticType')])
		names(disgenet.info) = c('pathway_id','pathway_name','pathway_class','pathway_semantic_type')
		
		hsap.genes = biomaRt::getBM(
			attributes=c('ensembl_gene_id','entrezgene_id'),
			mart = hsap)
		
		disgenet.all = disgenet.all[,c('geneId','diseaseId','diseaseClass','score','source')]
		names(disgenet.all) = c('entrezgene_id','pathway_id','pathway_class','score','source')
		
		disgenet.all = merge(disgenet.all,hsap.genes,by='entrezgene_id')
		
		if (!is.null(disgenet.class) & !identical(disgenet.class,'all')) {
			disgenet.classes = do.call(rbind,if (is.null(n.cores)) {
				lapply(1:nrow(disgenet.all),function(i) if (nchar(disgenet.all[i,]$pathway_class)) data.frame(disgenet.all[i,][setdiff(colnames(disgenet.all),'pathway_class')],pathway_class=unlist(strsplit(disgenet.all[i,]$pathway_class,';'))) else data.frame(disgenet.all[i,][setdiff(colnames(disgenet.all),'pathway_class')],pathway_class=''))
			} else {
				parallel::mclapply(1:nrow(disgenet.all),function(i) if (nchar(disgenet.all[i,]$pathway_class)) data.frame(disgenet.all[i,][setdiff(colnames(disgenet.all),'pathway_class')],pathway_class=unlist(strsplit(disgenet.all[i,]$pathway_class,';'))) else data.frame(disgenet.all[i,][setdiff(colnames(disgenet.all),'pathway_class')],pathway_class=''),mc.cores=n.cores)
			})
			n.before = nrow(unique(disgenet.all[c('pathway_id','ensembl_gene_id')]))
			disgenet.all = unique(subset(disgenet.classes,pathway_class %in% disgenet.class)[setdiff(colnames(disgenet.classes),'pathway_class')])
			n.after = nrow(unique(disgenet.all[c('pathway_id','ensembl_gene_id')]))
			message(paste0('Removed ',n.before-n.after,' annotations not in class(es): "',paste(disgenet.class,collapse='", "'),'"'))
			
			database = paste0(database,' (',paste(disgenet.class,collapse=', '),')')
		}
		if (!is.null(disgenet.source) & !identical(disgenet.source,'all')) {
			disgenet.sources = do.call(rbind,if (is.null(n.cores)) {
				lapply(1:nrow(disgenet.all),function(i) if (nchar(disgenet.all[i,]$source)) data.frame(disgenet.all[i,][setdiff(colnames(disgenet.all),'source')],source=unlist(strsplit(disgenet.all[i,]$source,';'))) else data.frame(disgenet.all[i,][setdiff(colnames(disgenet.all),'source')],source=''))
			} else {
				parallel::mclapply(1:nrow(disgenet.all),function(i) if (nchar(disgenet.all[i,]$source)) data.frame(disgenet.all[i,][setdiff(colnames(disgenet.all),'source')],source=unlist(strsplit(disgenet.all[i,]$source,';'))) else data.frame(disgenet.all[i,][setdiff(colnames(disgenet.all),'source')],source=''),mc.cores=n.cores)
			})
			n.before = nrow(unique(disgenet.all[c('pathway_id','ensembl_gene_id')]))
			disgenet.all = unique(subset(disgenet.sources,source %in% disgenet.source)[setdiff(colnames(disgenet.sources),'source')])
			n.after = nrow(unique(disgenet.all[c('pathway_id','ensembl_gene_id')]))
			message(paste0('Removed ',n.before-n.after,' annotations not from source(s): "',paste(disgenet.source,collapse='", "'),'"'))
			
			database = paste0(database,' (',paste(disgenet.source,collapse=', '),')')
		}
		if (!is.null(disgenet.score)) {
			n.before = nrow(unique(disgenet.all[c('pathway_id','ensembl_gene_id')]))
			disgenet.all = subset(disgenet.all,score >= disgenet.score)
			n.after = nrow(unique(disgenet.all[c('pathway_id','ensembl_gene_id')]))
			message(paste0('Removed ',n.before-n.after,' annotations with score < ',disgenet.score))
		}
		disgenet.all = unique(disgenet.all[,c('ensembl_gene_id','pathway_id')])
		
		gene2disgenet = if (is.null(n.cores)) {
			lapply(split(disgenet.all,disgenet.all$ensembl_gene_id),function(x) unique(x$pathway_id))
		} else {
			parallel::mclapply(split(disgenet.all,disgenet.all$ensembl_gene_id),function(x) unique(x$pathway_id),mc.cores=n.cores)
		}
		
		disgenet.info = subset(disgenet.info,pathway_id %in% unique(disgenet.all$pathway_id))
		
		rownames(disgenet.info) = NULL
		result = createPathwayAnnotationObject(
			database = database,
			organism = ensembl.species,
			annotation = gene2disgenet,
			definition = disgenet.info,
			ensembl.version = ensembl.version,
			annotation.version = disgenet.version
		)
	} else if (toupper(database) == 'REACTOMEPPI') {
		database = 'ReactomePPI'
		
		temp.file = tempfile() # temp.file='human_annotated_PPIs.txt.gz'

		if (reactome.species == 'HSA') {
			download.file('https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.tab-delimited.txt',temp.file,quiet=TRUE,timeout=600)
		} else {
			download.file('https://reactome.org/download/current/interactors/reactome.all_species.interactions.tab-delimited.txt',temp.file,quiet=TRUE,timeout=600)
		}
		reactomeppi.all = read.delim(
			temp.file,
			header=FALSE,
			skip=1,
			col.names=c('pathway1','ensembl1','entrez1','pathway2','ensembl2','entrez2','interaction_type','interaction_context','pubmed'))
		reactomeppi.pathway.version = paste('Reactome',rbioapi::rba_reactome_version(verbose=FALSE))
		
		reactomeppi.original = reactomeppi.all
		
		reactomeppi.columns = c('pathway1','ensembl1','pathway2','ensembl2','interaction_type','interaction_context','pubmed')
		reactomeppi.all = subset(reactomeppi.all,select=reactomeppi.columns)
		
		reactomeppi.all$interaction_id = unlist(lapply(1:nrow(reactomeppi.all),function(i) {
			paste(sort(c(reactomeppi.all[i,'pathway1'],reactomeppi.all[i,'pathway2'])),collapse='|')
		}))
		
		# Make mirror image, switching proteins in columns 1 and 2 (and flipping directionality in the meantime)
		reactomeppi.dir1 = c('pathway1','ensembl2','interaction_type','interaction_context','pubmed')
		reactomeppi.dir2 = c('pathway2','ensembl1','interaction_type','interaction_context','pubmed')

		reactomeppi.duplicated = subset(reactomeppi.all,pathway1 == pathway2 & !duplicated(interaction_id),select=reactomeppi.dir1)
		reactomeppi.unique = subset(reactomeppi.all,pathway1 != pathway2 & !duplicated(interaction_id))
		
		reactomeppi.1 = subset(reactomeppi.unique,select=reactomeppi.dir1)
		reactomeppi.2 = subset(reactomeppi.unique,select=reactomeppi.dir2)
		
		names(reactomeppi.1) = names(reactomeppi.2) = names(reactomeppi.duplicated) = c('pathway_id','ensembl_id','interaction_type','interaction_context','pubmed')

		reactomeppi.all = rbind(reactomeppi.duplicated,reactomeppi.1,reactomeppi.2)
	
		gene2reactomeppi.df = if (is.null(n.cores)) {
			do.call(rbind,lapply(1:nrow(reactomeppi.all),function(i) {
				x = reactomeppi.all[i,]
				y = data.frame(
					pathway_id = x$pathway_id,
					ensembl_id = gsub('^ENSEMBL:','',unlist(strsplit(x$ensembl_id,'\\|'))),
					x[c('interaction_type','interaction_context','pubmed')]
				)
				y
			}))
		} else {
			do.call(rbind,parallel::mclapply(1:nrow(reactomeppi.all),function(i) {
				x = reactomeppi.all[i,]
				y = data.frame(
					pathway_id = x$pathway_id,
					ensembl_id = gsub('^ENSEMBL:','',unlist(strsplit(x$ensembl_id,'\\|'))),
					x[c('interaction_type','interaction_context','pubmed')]
				)
				y
			},mc.cores=n.cores))
		}
		
		gene2reactomeppi.df$ensembl_id = gsub('\\.[0-9]+$','',gene2reactomeppi.df$ensembl_id)
		
		gene2reactomeppi.df = subset(
			gene2reactomeppi.df,
			grepl('^ENSG[0-9]{11}$|^ENST[0-9]{11}$|^ENSP[0-9]{11}$',ensembl_id)
		)
		
		ens.genes = biomaRt::getBM(
			attributes=c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id'),
			mart = ens)
		
		ens.transcripts = subset(ens.genes,grepl('^ENST[0-9]{11}$',ensembl_transcript_id))
		ens.peptides = subset(ens.genes,grepl('^ENSP[0-9]{11}$',ensembl_peptide_id))
		
		enst.to.ensg = ens.transcripts$ensembl_gene_id
		names(enst.to.ensg) = ens.transcripts$ensembl_transcript_id
		ensp.to.ensg = ens.peptides$ensembl_gene_id
		names(ensp.to.ensg) = ens.peptides$ensembl_peptide_id
		
		gene2reactomeppi.df$ensembl_gene_id = with(gene2reactomeppi.df,ifelse(
			grepl('^ENSG[0-9]{11}',ensembl_id),
			ensembl_id,
			ifelse(
				grepl('^ENST[0-9]{11}',ensembl_id),
				enst.to.ensg[ensembl_id],
				ifelse(
					grepl('^ENSP[0-9]{11}',ensembl_id),
					ensp.to.ensg[ensembl_id],
					''
				)
			)
		))
		gene2reactomeppi.df = unique(subset(gene2reactomeppi.df,select=c('pathway_id','ensembl_gene_id','interaction_type','interaction_context','pubmed')))
				
		reactomeppi.info = unique(gene2reactomeppi.df['pathway_id'])
		# table(gsub(':.+','',reactomeppi.info$pathway_id))
		
		# uniprotkb
		uniprot.names = UniProt.ws::mapUniProt(
			from = 'UniProtKB_AC-ID',
			to = 'Gene_Name',
			query = unique(gsub('uniprotkb:','',subset(reactomeppi.info,grepl('uniprotkb',pathway_id))$pathway_id)),
			verbose = FALSE
		)
		
		uniprot.names = do.call(rbind,lapply(split(uniprot.names, uniprot.names$From),function(x) data.frame(uniprot_id = unique(x$From),name=paste(x$To,collapse='|'))))
		
		uniprot.to.name = uniprot.names$name
		names(uniprot.to.name) = paste0('uniprotkb:',uniprot.names$uniprot_id)
		
		# reactome	
		reactome.1 = subset(reactomeppi.original,grepl('^reactome:',pathway1),select=c('pathway1','ensembl1'))
		reactome.2 = subset(reactomeppi.original,grepl('^reactome:',pathway2),select=c('pathway2','ensembl2'))
		names(reactome.1) = names(reactome.2) = c('reactome_gene_id','ensembl')
		reactome.names = subset(unique(rbind(reactome.1,reactome.2)),reactome_gene_id %in% reactomeppi.info$pathway_id)
		
		reactome.names$ensembl = gsub('ENSEMBL:','',reactome.names$ensembl)
		reactome.names = subset(reactome.names,nchar(ensembl) > 1)

		ens.gene.names = biomaRt::getBM(
			attributes=c('ensembl_gene_id','external_gene_name'),
			mart = ens)
		genes.to.symbol = ens.gene.names$external_gene_name
		names(genes.to.symbol) = ens.gene.names$ensembl_gene_id
		
		reactome.names$reactome_gene_name = genes.to.symbol[reactome.names$ensembl]
		reactome.names = subset(reactome.names,!is.na(reactome.names$reactome_gene_name))
		
		reactome.to.name = reactome.names$reactome_gene_name
		names(reactome.to.name) = reactome.names$reactome_gene_id
		
		# ChEBI
		temp.file = tempfile()
		download.file('https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz',temp.file,quiet=TRUE,timeout=600)
		file.rename(temp.file,paste0(temp.file,'txt.gz'))
		temp.file=paste0(temp.file,'txt.gz')
		
		chebi.compounds = data.table::fread(temp.file,data.table=FALSE)
		names(chebi.compounds) = c('chebi_id','status','chebi_compound_id','source','parent_id','name','definition','modified','created','star')
		
		# chebi.compounds$chebi_compound_id = paste0('ChEBI:',format(chebi.compounds$chebi_compound_id,scientific=FALSE,trim=TRUE))
		chebi.compounds$chebi_compound_id = gsub('CHEBI','ChEBI',chebi.compounds$chebi_compound_id)
		
		chebi.to.name = chebi.compounds$name
		names(chebi.to.name) = chebi.compounds$chebi_compound_id
		
		chebi.to.description = gsub('^null$','',chebi.compounds$definition)
		names(chebi.to.description) = chebi.compounds$chebi_compound_id
		
		# reactomeppi.info = subset(reactomeppi.info,
		# 	pathway_id %in% c(names(uniprot.to.name),names(reactome.to.name),names(chebi.to.name))
		# )
		
		reactomeppi.info$pathway_name = with(reactomeppi.info,ifelse(
			grepl('^uniprotkb:',pathway_id),
			uniprot.to.name[pathway_id],
			ifelse(
				grepl('^reactome:',pathway_id),
				reactome.to.name[pathway_id],
				ifelse(
					grepl('^ChEBI:',pathway_id),
					chebi.to.name[pathway_id],
					''
				)
			)
		))
		reactomeppi.info$pathway_name[is.na(reactomeppi.info$pathway_name)] = reactomeppi.info$pathway_id[is.na(reactomeppi.info$pathway_name)]
		reactomeppi.info$pathway_definition = with(reactomeppi.info,ifelse(
			grepl('^ChEBI:',pathway_id),
				chebi.to.description[pathway_id],
				''
			)
		)
				
		gene2reactomeppi = if (is.null(n.cores)) {
			lapply(split(gene2reactomeppi.df,gene2reactomeppi.df$ensembl_gene_id),function(x) unique(x$pathway_id))
		} else {
			parallel::mclapply(split(gene2reactomeppi.df,gene2reactomeppi.df$ensembl_gene_id),function(x) unique(x$pathway_id),mc.cores=n.cores)
		}
		
		reactomeppi.info = subset(reactomeppi.info,pathway_id %in% unique(reactomeppi.all$pathway_id))
		
		rownames(reactomeppi.info) = NULL
		result = createPathwayAnnotationObject(
			database = database,
			organism = ensembl.species,
			annotation = gene2reactomeppi,
			definition = reactomeppi.info,
			ensembl.version = ensembl.version,
			annotation.version = reactomeppi.pathway.version
		)
	} else if (toupper(database) %in% c('IID','PPI','OPHID')) {
		database = 'IID'
		
		temp.file = tempfile() # temp.file='human_annotated_PPIs.txt.gz'
		# download.file(paste0('http://ophid.utoronto.ca/static/download/',iid.species,'_annotated_PPIs.txt.gz'),temp.file,quiet=TRUE)
		download.file(paste0('http://iid.ophid.utoronto.ca/static/download/',iid.species,'_annotated_PPIs.txt.gz'),temp.file,quiet=TRUE,timeout=600)
		iid.all = read.delim(
			gzfile(temp.file),
			header=TRUE)
		iid.version = paste0('Downloaded ',`attr<-`(as.POSIXct(Sys.time()),'tzone','UTC'),' UTC')
		
		iid.columns = c('uniprot1','uniprot2','symbol1','symbol2','methods','pmids','db_with_ppi','evidence_type','n_methods','n_exp_pmids','n_pred_pmids','n_pmids','evidence_exp','evidence_ortho','evidence_comp','directed','bidirected','directions','stable','transient')
		iid.all = subset(iid.all,select=iid.columns)
		
		# Make mirror image, switching proteins in columns 1 and 2 (and flipping directionality in the meantime)
		iid.flipped = c('uniprot2','uniprot1','symbol2','symbol1','methods','pmids','db_with_ppi','evidence_type','n_methods','n_exp_pmids','n_pred_pmids','n_pmids','evidence_exp','evidence_ortho','evidence_comp','directed','bidirected','directions','stable','transient')
		func = function(x) {
			this.direction = unique(x$directions)
			x = x[,iid.flipped] # Rearrange columns
			names(x) = iid.columns
			if (this.direction %in% c('>','<')) { # If the PPI was directional, flip the direction
				x$directions = with(x,ifelse(uniprot1==uniprot2,'><',ifelse(this.direction=='>','<','>'))) # If evidence was directional but protein interacts with self, mark as bidirectional instead
			}
			x
		}
		iid.mirror = do.call(rbind,if (is.null(n.cores)) {
			lapply(split(iid.all,iid.all$directions),func)
		} else {
			parallel::mclapply(split(iid.all,iid.all$directions),func,mc.cores=n.cores)
		})
		iid.all = rbind(subset(iid.all,uniprot1 != uniprot2),iid.mirror)
		
		ens.genes = biomaRt::getBM(
			attributes=c('ensembl_gene_id','uniprotswissprot','description'),
			mart = ens)
		
		iid.all = merge(iid.all,ens.genes,by.x='uniprot1',by.y='uniprotswissprot',all.x=TRUE)
		
		iid.info = unique(iid.all[,c('uniprot1','symbol1','description')])
		names(iid.info) = c('pathway_id','pathway_name','pathway_description')
		
		iid.info$pathway_description[is.na(iid.info$pathway_description)] = ''
		
		if (any(duplicated(iid.info$pathway_id))) {
			duplicated.pathways = iid.info[duplicated(iid.info$pathway_id),'pathway_id']
			iid.duplicated = subset(iid.info,pathway_id %in% duplicated.pathways)
			iid.duplicated = do.call(rbind,if (is.null(n.cores)) {
				lapply(split(iid.duplicated,iid.duplicated$pathway_id),function(x) data.frame(unique(x[,c('pathway_id','pathway_name')]),pathway_description=paste(x$pathway_description,collapse='|')))
			} else {
				parallel::mclapply(split(iid.duplicated,iid.duplicated$pathway_id),function(x) data.frame(unique(x[,c('pathway_id','pathway_name')]),pathway_description=paste(x$pathway_description,collapse='|')),mc.cores=n.cores)
			})
			iid.info = rbind(subset(iid.info,!pathway_id %in% duplicated.pathways),iid.duplicated)
			iid.info = iid.info[order(iid.info$pathway_id),]
			rownames(iid.info) = NULL
		}

		iid.all = subset(iid.all,!is.na(ensembl_gene_id))[,c('ensembl_gene_id','uniprot2','methods','db_with_ppi','evidence_type','n_methods','n_pmids','directions','stable','transient')]
		names(iid.all)[names(iid.all) == 'uniprot2'] = 'pathway_id'

		if (!is.null(iid.n.methods)) {
			n.before = nrow(unique(iid.all[c('pathway_id','ensembl_gene_id')]))
			iid.all = subset(iid.all,n_methods >= iid.n.methods)
			n.after = nrow(unique(iid.all[c('pathway_id','ensembl_gene_id')]))
			message(paste0('Removed ',n.before-n.after,' annotations with number of methods < ',iid.n.methods))
		}

		if (!is.null(iid.n.publications)) {
			n.before = nrow(unique(iid.all[c('pathway_id','ensembl_gene_id')]))
			iid.all = subset(iid.all,n_pmids >= iid.n.publications)
			n.after = nrow(unique(iid.all[c('pathway_id','ensembl_gene_id')]))
			message(paste0('Removed ',n.before-n.after,' annotations with number of publications (PMID) < ',iid.n.publications))
		}

		if (!is.null(iid.evidence.type) & !identical(iid.evidence.type,'all')) {
			iid.evidence.types = do.call(rbind,if (is.null(n.cores)) {
				lapply(1:nrow(iid.all),function(i) if (nchar(iid.all[i,]$evidence_type)) data.frame(iid.all[i,][setdiff(colnames(iid.all),'evidence_type')],evidence_type=unlist(strsplit(iid.all[i,]$evidence_type,'\\|'))) else data.frame(iid.all[i,][setdiff(colnames(iid.all),'evidence_type')],evidence_type=''))
			} else {
				parallel::mclapply(1:nrow(iid.all),function(i) if (nchar(iid.all[i,]$evidence_type)) data.frame(iid.all[i,][setdiff(colnames(iid.all),'evidence_type')],evidence_type=unlist(strsplit(iid.all[i,]$evidence_type,'\\|'))) else data.frame(iid.all[i,][setdiff(colnames(iid.all),'evidence_type')],evidence.type=''),mc.cores=n.cores)
			})
			n.before = nrow(unique(iid.all[c('pathway_id','ensembl_gene_id')]))
			iid.all = unique(subset(iid.evidence.types,evidence_type %in% iid.evidence.type)[setdiff(colnames(iid.all),'evidence_type')])
			n.after = nrow(unique(iid.all[c('pathway_id','ensembl_gene_id')]))
			message(paste0('Removed ',n.before-n.after,' annotations not supported by evidence type(s): "',paste(iid.evidence.type,collapse='", "'),'"'))
			
			database = paste0(database,' (',paste(iid.evidence.type,collapse=', '),')')
		}

		if (!is.null(iid.direction) & !identical(iid.direction,'all')) {
			if (length(iid.direction) > 1) {
				iid.direction = iid.direction[1]
				warning('iid.direction has length > 1. Using only the first value.')
			}
			iid.direction = if (iid.direction %in% c('>','<')) {
				direction.text = if (iid.direction == '>') 'directional gene->protein' else if (iid.direction == '<') 'directional gene<-protein'
				iid.direction = c(iid.direction,'><')
			} else if (iid.direction == '><') {
				direction.text = 'bidirectional gene<->protein'
			}
			n.before = nrow(unique(iid.all[c('pathway_id','ensembl_gene_id')]))
			iid.all = subset(iid.all,directions %in% iid.direction)[setdiff(colnames(iid.all),'directions')]
			n.after = nrow(unique(iid.all[c('pathway_id','ensembl_gene_id')]))
			message(paste0('Removed ',n.before-n.after,' annotations without the following type(s) of directional evidence: "',paste(iid.direction,collapse='", "'),'"'))
			
			database = paste0(database,' (',paste(direction.text,collapse=', '),')')
		}
		
		iid.all = unique(iid.all[,c('ensembl_gene_id','pathway_id')])
		
		gene2iid = if (is.null(n.cores)) {
			lapply(split(iid.all,iid.all$ensembl_gene_id),function(x) unique(x$pathway_id))
		} else {
			parallel::mclapply(split(iid.all,iid.all$ensembl_gene_id),function(x) unique(x$pathway_id),mc.cores=n.cores)
		}
		
		iid.info = subset(iid.info,pathway_id %in% unique(iid.all$pathway_id))
		
		rownames(iid.info) = NULL
		result = createPathwayAnnotationObject(
			database = database,
			organism = ensembl.species,
			annotation = gene2iid,
			definition = iid.info,
			ensembl.version = ensembl.version,
			annotation.version = iid.version
		)
	}
	result
}

#' Lift an annotation over to another organism on the basis of orthology information.
#'
#' This function provides a workflow for converting an annotation for use with a new
#' organism. An internet connection is required
#' 
#' @param annotation An object of pathwayAnnotation class (query species).
#' @param ensembl.species A character string specifying the organism for which genes (ENSEMBL) should be annotated (target species). Valid values can be found using `fetchOrganisms("ENSEMBL")`.
#' @param ensembl.version A character string specifying the version of the ENSEMBL database from which to derive annotations for the target species.
#' @param orthology.type A character string specifying the type(s) of orthology allowed for transferring annotations to the target species. Possible values are `"ortholog_one2one"`, `"ortholog_one2many"`, `"ortholog_many2many"`. If `NULL` (default), `"ortholog_one2one"` (only single-copy orthologs) is assumed.
#' @param percent.identity.1 A number specifying the percent identity threshold (target to query) for orthologs to include. Valid values are from `0` to `100`. If `NULL`, no filter is applied.
#' @param percent.identity.2 A number specifying the percent identity threshold (query to target) for orthologs to include. Valid values are from `0` to `100`. If `NULL`, no filter is applied.
#' @param goc.score A number specifying the gene order conservation score threshold for orthologs to include. Valid values are from `0` to `100`. If `NULL`, no filter is applied.
#' @param wga.coverage A number specifying the whole genome alignment score threshold for orthologs to include. Valid values are from `0` to `100`. If `NULL`, no filter is applied.
#' @param orthology.confidence A number specifying the orthology confidence threshold for orthologs to include. Valid values are from `0` to `1`. If `NULL`, no filter is applied.
#' @param n.cores A number specifying the number of cores to use.
#' @return An object of pathwayAnnotation class
#' @export
liftAnnotations = function(
	annotation,
	ensembl.species,
	ensembl.version = NULL,
	orthology.type = NULL,
	percent.identity.1 = NULL, # % identity target to query
	percent.identity.2 = NULL, # % identity query to target,
	goc.score = NULL,
	wga.coverage = NULL,
	orthology.confidence = NULL,
	n.cores = NULL
) {
	# suppressPackageStartupMessages(require(biomaRt))
	# if (!is.null(n.cores)) suppressPackageStartupMessages(require(parallel))

	from.species = annotation@organism
	to.species = ensembl.species
	if (from.species == to.species) stop('Target and query species cannot be identical.')
	
	if (is.null(ensembl.version)) {
		ens = biomaRt::useEnsembl(
			biomart = 'ENSEMBL_MART_ENSEMBL',
			dataset=paste0(from.species,'_gene_ensembl')
		)
		ensembl.version = gsub('Ensembl Genes (.+)','Ensembl \\1',subset(biomaRt::listEnsembl(),biomart == 'genes')$version)
	} else {
		ensembl.version = as.character(ensembl.version)
		ens = biomaRt::useEnsembl(
			biomart = 'ENSEMBL_MART_ENSEMBL',
			version = ensembl.version,
			dataset=paste0(from.species,'_gene_ensembl')
		)
		ensembl.version = subset(biomaRt::listEnsemblArchives(),version == ensembl.version)$name
	}
	
	ens.attr = biomaRt::listAttributes(ens)
	good.attr = intersect(ens.attr$name,c(
			'ensembl_gene_id',
			paste0(to.species,'_homolog_ensembl_gene'),
			paste0(to.species,'_homolog_orthology_type'),
			paste0(to.species,'_homolog_perc_id'),
			paste0(to.species,'_homolog_perc_id_r1'),
			paste0(to.species,'_homolog_goc_score'),
			paste0(to.species,'_homolog_wga_coverage'),
			paste0(to.species,'_homolog_orthology_confidence')
		)
	)
	ens.genes = biomaRt::getBM(
		attributes=good.attr,
		mart = ens)
	
	ens.genes = ens.genes[!is.na(ens.genes[[paste0(to.species,'_homolog_ensembl_gene')]]) & as.logical(nchar(ens.genes[[paste0(to.species,'_homolog_ensembl_gene')]])),]
	
	ens.genes = subset(ens.genes,ensembl_gene_id %in% names(annotation@annotation))
	
	message('Matched ',length(na.omit(unique(ens.genes[[paste0(to.species,'_homolog_ensembl_gene')]]))),' target genes.')
	
	# Filter
	if (is.null(orthology.type)) {
		filter.genes = ens.genes[[paste0(to.species,'_homolog_orthology_type')]] %in% 'ortholog_one2one'
		message(paste0('Removing ',length(na.omit(unique(ens.genes[!filter.genes,paste0(to.species,'_homolog_ensembl_gene')]))),' target genes via orthology filter. Keeping ortholog type(s): "ortholog_one2one".'))
		ens.genes = ens.genes[filter.genes,]
	} else if (orthology.type == 'all') {
		message(paste0('Removing 0 target genes via orthology filter. Keeping ortholog type(s): "',paste(unique(ens.genes[[paste0(to.species,'_homolog_orthology_type')]]),collapse='", "'),'".'))
	} else if (length(orthology.type) && class(orthology.type) == 'character') {
		filter.genes = ens.genes[[paste0(to.species,'_homolog_orthology_type')]] %in% orthology.type
		message(paste0('Removing ',length(na.omit(unique(ens.genes[!filter.genes,paste0(to.species,'_homolog_ensembl_gene')]))),' target genes via orthology filter. Keeping ortholog type(s): "',paste(orthology.type,collapse='", "'),'".'))
		ens.genes = ens.genes[filter.genes,]
	} else {
		stop('Error')
	}
	if (!is.null(percent.identity.1)) {
		if (paste0(to.species,'_homolog_perc_id') %in% good.attr) {
			filter.genes = !is.na(ens.genes[[paste0(to.species,'_homolog_perc_id')]]) & ens.genes[[paste0(to.species,'_homolog_perc_id')]] >= percent.identity.1
			message(paste0('Removing ',length(na.omit(unique(ens.genes[!filter.genes,paste0(to.species,'_homolog_ensembl_gene')]))),' target genes where percent identity < ',percent.identity.1,'% (',to.species,' query to ',from.species,' target).'))
			ens.genes = ens.genes[filter.genes,]
		} else {
			warning(paste0('Not filtering on "percent.identity.1" as attribute "',paste0(to.species,'_homolog_perc_id'),'" is not available via bioMart.'))
		}
	}
	if (!is.null(percent.identity.2)) {
		if (paste0(to.species,'_homolog_perc_id_r1') %in% good.attr) {
			filter.genes = !is.na(ens.genes[[paste0(to.species,'_homolog_perc_id_r1')]]) & ens.genes[[paste0(to.species,'_homolog_perc_id_r1')]] >= percent.identity.2
			message(paste0('Removing ',length(na.omit(unique(ens.genes[!filter.genes,paste0(to.species,'_homolog_ensembl_gene')]))),' target genes where percent identity < ',percent.identity.2,'% (',from.species,' target to ',to.species,' query).'))
			ens.genes = ens.genes[filter.genes,]
		} else {
			warning(paste0('Not filtering on "percent.identity.2" as attribute "',paste0(to.species,'_homolog_perc_id_r1'),'" is not available via bioMart.'))
		}
	}
	if (!is.null(goc.score)) {
		if (paste0(to.species,'_homolog_goc_score') %in% good.attr) {
			filter.genes = !is.na(ens.genes[[paste0(to.species,'_homolog_goc_score')]]) & ens.genes[[paste0(to.species,'_homolog_goc_score')]] >= goc.score
			message(paste0('Removing ',length(na.omit(unique(ens.genes[!filter.genes,paste0(to.species,'_homolog_ensembl_gene')]))),' target genes where gene-order conservation score < ',goc.score,'%.'))
			ens.genes = ens.genes[filter.genes,]
		} else {
			warning(paste0('Not filtering on "goc.score" as attribute "',paste0(to.species,'_homolog_goc_score'),'" is not available via bioMart.'))
		}
	}
	if (!is.null(wga.coverage)) {
		if (paste0(to.species,'_homolog_wga_coverage') %in% good.attr) {
			filter.genes = !is.na(ens.genes[[paste0(to.species,'_homolog_wga_coverage')]]) & ens.genes[[paste0(to.species,'_homolog_wga_coverage')]] >= wga.coverage
			message(paste0('Removing ',length(na.omit(unique(ens.genes[!filter.genes,paste0(to.species,'_homolog_ensembl_gene')]))),' target genes where whole-genome alignment coverage < ',wga.coverage,'%.'))
			ens.genes = ens.genes[filter.genes,]
		} else {
			warning(paste0('Not filtering on "wga.coverage" as attribute "',paste0(to.species,'_homolog_wga_coverage'),'" is not available via bioMart.'))
		}
	}
	if (!is.null(orthology.confidence)) {
		if (paste0(to.species,'_homolog_orthology_confidence') %in% good.attr) {
			filter.genes = !is.na(ens.genes[[paste0(to.species,'_homolog_orthology_confidence')]]) & ens.genes[[paste0(to.species,'_homolog_orthology_confidence')]] >= orthology.confidence
			message(paste0('Removing ',length(na.omit(unique(ens.genes[!filter.genes,paste0(to.species,'_homolog_ensembl_gene')]))),' target genes where orthology confidence < ',orthology.confidence,'.'))
			ens.genes = ens.genes[filter.genes,]
		} else {
			warning(paste0('Not filtering on "orthology.confidence" as attribute "',paste0(to.species,'_homolog_orthology_confidence'),'" is not available via bioMart.'))
		}
	}
	message('Matched ',length(unique(ens.genes[[paste0(to.species,'_homolog_ensembl_gene')]])),' target genes after filters.')
	
	pathway.all = if (is.null(n.cores)) {
		do.call(rbind,lapply(names(annotation@annotation),function(x) data.frame(ensembl_gene_id=x,pathway_id=annotation@annotation[[x]])))
	} else {
		do.call(rbind,parallel::mclapply(names(annotation@annotation),function(x) data.frame(ensembl_gene_id=x,pathway_id=annotation@annotation[[x]]),mc.cores=n.cores))
	}
	
	# Update gene:pathway associations
	pathway.all = unique(merge(pathway.all,ens.genes[,c('ensembl_gene_id',paste0(to.species,'_homolog_ensembl_gene'))],by='ensembl_gene_id')[,c(paste0(to.species,'_homolog_ensembl_gene'),'pathway_id')])
	names(pathway.all) = c('ensembl_gene_id','pathway_id')

	gene2pathway = if (is.null(n.cores)) {
		lapply(split(pathway.all,pathway.all$ensembl_gene_id),function(x) unique(x$pathway_id))
	} else {
		parallel::mclapply(split(pathway.all,pathway.all$ensembl_gene_id),function(x) unique(x$pathway_id),mc.cores=n.cores)
	}

	pathway.info = subset(annotation@definition,pathway_id %in% unique(pathway.all$pathway_id))
	
	result = createPathwayAnnotationObject(
		database = annotation@source,
		organism = to.species,
		annotation = gene2pathway,
		definition = pathway.info,
		ensembl.version = ensembl.version,
		annotation.version = if (grepl('\\(lifted from',annotation@pathway_annotation_version)) {
			paste0(gsub('\\((lifted from .+?)\\)$','(\\1, then ',annotation@pathway_annotation_version),from.species,')')
		} else {
			paste0(annotation@pathway_annotation_version,' (lifted from ',from.species,')')
		}
	)
}

#' Export an annotation into GMT (gene matrix transposed) file format.
#'
#' This function provides a workflow for exported an annotation as a GMT file. Exported files
#' may be used with tools such as g:Profiler.
#' 
#' @param annotation An object of pathwayAnnotation class.
#' @param description If character class, column (of definition object) from which to populate the optional description file. If `NULL`, the column will be populated with dummy ("na") values.
#' @param file File path for exported GMT file.
#' @param n.cores A number specifying the number of cores to use.
#' @export

exportGMT = function(
	annotation,
	description = 'pathway_name', # Column to populate the description of GMT file. Set to NULL if no description
	file = 'annotation_export.gmt',
	n.cores = NULL
) {
	if (!'pathwayAnnotation' %in% class(annotation)) stop('Annotation class must be "pathwayAnnotation"')
	
	annotation.list = getAnnotationsList(annotation)
	if (is.null(n.cores)) {
		annotation.df = do.call(rbind,lapply(names(annotation),function(x) {
			data.frame(ensembl_gene_id = x, pathway_id = annotation.list[[x]])
		}))
	} else {
		annotation.df = do.call(rbind,parallel::mclapply(names(annotation),function(x) {
			data.frame(ensembl_gene_id = x, pathway_id = annotation.list[[x]])
		},mc.cores=n.cores))
	}
	annotation.df = annotation.df[complete.cases(annotation.df),]
	
	if (is.null(description)) {
		pathway.df = getDefinitions(annotation)['pathway_id']
		pathway.df[['description']] = 'na'
	} else {
		pathway.df = getDefinitions(annotation)[c('pathway_id',description)]
		names(pathway.df)[2] = 'description'
	}
	pathway.df = subset(pathway.df,pathway_id %in% annotation.df$pathway_id)
	pathway.df = pathway.df[order(pathway.df$pathway_id),]
	rownames(pathway.df) = pathway.df$pathway_id
	
	annotation.split = split(annotation.df,annotation.df$pathway_id)
	if (is.null(n.cores)) {
		gmt.out = unlist(lapply(pathway.df$pathway_id,function(x) {
			out = c(unlist(pathway.df[x,]),annotation.split[[x]]$ensembl_gene_id)
			names(out) = NULL
			paste(out,collapse='\t')
		}))
	} else {
		gmt.out = unlist(parallel::mclapply(pathway.df$pathway_id,function(x) {
			out = c(unlist(pathway.df[x,]),annotation.split[[x]]$ensembl_gene_id)
			names(out) = NULL
			paste(out,collapse='\t')
		},mc.cores=n.cores))
	}
	write(gmt.out,file=file,sep='\n')
}

#' Import an annotation from GMT (gene matrix transposed) file format.
#'
#' This function provides a workflow for importing an annotation as a GMT file.
#' 
#' @param file File path for GMT file.
#' @param database Name of the annotation source.
#' @param organism Name of the organism. It is advisable to use ENSEMBL codes, which can be retrieved using `fetchOrganisms('ENSEMBL')`.
#' @param n.cores A number specifying the number of cores to use.
#' @param ... Additional arguments to be passed to `createPathwayAnnotationObject`.
#' @export

importGMT = function(
	file,
	database,
	organism,
	n.cores = NULL,
	...
) {
	gmt.in = strsplit(scan(file, what='', sep='\n', quiet=TRUE),'\t')
	if (is.null(n.cores)) {
		gmt.df = do.call(rbind,lapply(gmt.in,function(x) {
			data.frame(pathway_id = x[1], pathway_name = x[2], ensembl_gene_id = x[3:length(x)])
		}))
	} else {
		gmt.df = do.call(rbind,parallel::mclapply(gmt.in,function(x) {
			data.frame(pathway_id = x[1], pathway_name = x[2], ensembl_gene_id = x[3:length(x)])
		},mc.cores=n.cores))
	}
	definition = unique(gmt.df[c('pathway_id','pathway_name')])
	
	gene2pathway = split(gmt.df$pathway_id,gmt.df$ensembl_gene_id)
	
	createPathwayAnnotationObject(
		database = database,
		annotation = gene2pathway,
		definition = definition,
		organism = organism,
		...)
}
