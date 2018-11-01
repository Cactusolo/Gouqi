# Gouqi

## we are going to launch a new project insighted by [Wang et al., (2018)](http://rspb.royalsocietypublishing.org/content/285/1890/20181742?from=groupmessage).

## Task List
- [x] **prepare species list**
	* all the recorded names
	* validated by The Plant List (TPL) using [Taxonstand](https://cran.r-project.org/package=Taxonstand)
	@Ruuyu please check file "data/OpenTree_data/species_list/Gouqi_checkname_TPLQery.csv" and "data/Lycium_PL_no_name_TPLQery.csv"
	* **Outgroup**
		_maybe we need to select a few species from its closed relatives or tribes inside Solanoideae (e.g, Grabowskia, Nicandreae?)_
- [] **data mining**
1. **_morphological data_**
	* trait characters
	* Medical traits
2. **_molecular data_**
	* GenBank data mining using [PhylotaR](https://github.com/ropensci/phylotaR)
	* choose gene with most coverage
	* sequence length
3. **_distribution data_**
	* GBIF
	* iDigBio
	* NSII (CVH?)
4. **_Environmental data_**
	* [WorldClim](http://www.worldclim.org/)?

- [] **data cleaning and evaluation**
- [] **phylogeny and dating**
- [] **niche modelling _species distribution model_**
- [] **wrapping up with a paper**
		@Ruuyu, @Christmas-Wu, @Cactusolo

## other stuff
* lists and tables better use ".csv" format put into "data" folder
* references and other Documents put into "reference" folder
* likewise scripts go to "script" foder

