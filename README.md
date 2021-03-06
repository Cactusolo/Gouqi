# Gouqi

## we are going to launch a new project insighted by [Wang et al., (2018)](http://rspb.royalsocietypublishing.org/content/285/1890/20181742?from=groupmessage).

## Task List
- [x] **prepare species list**
	* all the recorded names
	* validated by The Plant List (TPL) using [Taxonstand](https://cran.r-project.org/package=Taxonstand)
	@Ruuyu please check file "data/OpenTree_data/species_list/Gouqi_checkname_TPLQery.csv" and "data/Lycium_PL_no_name_TPLQery.csv"
	* **Outgroup**
		Mybe we need to select a few species from its closed relatives or tribes inside Solanoideae (e.g, Grabowskia, Nicandreae?)
		+ _Atropa belladonna_
		+ _Jaborosa integrifolia_
		+ _Nolana arenicola_
	
- [] **data mining**
1. **_morphological data_**
	* trait characters
	* Medical traits
2. [x] **_molecular data_**
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
	 * adding constraint tree from [Miller et al. (2011)](https://academic.oup.com/mbe/article/28/1/793/986705) [x]
	 * updated the phylogeny with Phphlawd data and constraint tree from Miller et al. (2011) and the file path: \results\tree\updated_tree12032018
	 * two names were validated by the plant list: 
	 + _Lycium arenicola_ is the syn of _Lycium cinereum_
	 + _Lycium tenue_ is the syn of _Lycium acutifolium_

- [] **niche modelling _species distribution model_**
- [] **wrapping up with a paper**
		@Ruuyu, @Christmas-Wu, @Cactusolo

## other stuff
* lists and tables better use ".csv" format put into "data" folder
* references and other Documents put into "reference" folder
* likewise scripts go to "script" foder

