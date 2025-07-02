import numpy as np
import pandas as pd

""" GCAM land type groups and aggregations. """
forest_names = ['Forest', 'ProtectedUnmanagedForest', 'UnmanagedForest']
pasture_names = ['Pasture', 'UnmanagedPasture', 'ProtectedUnmanagedPasture']
grassland_names = ['Grassland', 'ProtectedGrassland']
crop_names_nonstandard = ['Corn', 'CornC4', 'FiberCrop', 'FodderGrass', 'FodderHerb', 'FodderHerbC4', 'Fruits', 'FruitsTree', 'Legumes', 'MiscCrop', 
	'MiscCropC4', 'MiscCropTree', 'NutsSeeds', 'NutsSeedsTree', 'OilCrop', 'OilCropTree', 'OilPalmTree', 'OtherGrain', 'OtherGrainC4', 'Rice', 
    'RootTuber', 'Soybean', 'SugarCrop', 'SugarCropC4', 'Vegetables', 'Wheat', 'biomass', 'biomassGrass', 'biomassTree'] # Nonstandard crop names.
crop_names = ['BioenergyCrop', 'Corn', 'FiberCrop', 'FodderGrass', 'FodderHerb', 'FodderHerb', 'Fruits', 'Legumes', 'MiscCrop', 'NutsSeeds',  
    'OilCrop', 'OilPalm', 'OtherGrain', 'Rice', 'RootTuber', 'Soybean', 'SugarCrop', 'Vegetables', 'Wheat']   # Standardized set of crop names.
other_arable_names = ['OtherArableLand']		# This is considered a managed cropland type.
shrubland_names = ['Shrubland', 'ProtectedShrubland']
urban_names = ['UrbanLand']                     # These are constant.
other_names = ['RockIceDesert', 'Tundra']       # These are constant.

""" Dictionary of GCAM land type groups and aggregations with standard crop names. """
gcam_landtype_groups = {'forest': forest_names, 'pasture': pasture_names, 'grass': grassland_names, 'crop': crop_names,
        'other_arable': other_arable_names, 'shrub': shrubland_names, 'urban': urban_names, 'other': other_names}

""" Dictionary of GCAM land type groups and aggregations with nonstandard crop names. """
gcam_landtype_groups_nonstandard = gcam_landtype_groups
gcam_landtype_groups_nonstandard['crop'] = crop_names_nonstandard

""" Dictionary of GCAM crop name mappings (keys = old names, values = new names) to create a standardized set of crop names across all files. """
gcam_crop_names = {'biomass': 'BioenergyCrop', 'biomassGrass': 'BioenergyCrop', 'biomassTree': 'BioenergyCrop', 'CornC4': 'Corn', 
        'FodderHerbC4': 'FodderHerb', 'FruitsTree': 'Fruits', 'MiscCropTree': 'MiscCrop', 'MiscCropC4': 'MiscCrop', 'NutsSeedsTree': 'NutsSeeds',
        'OilCropTree': 'OilCrop', 'OilPalmTree': 'OilPalm', 'OtherGrainC4': 'OtherGrain', 'SugarCropC4': 'SugarCrop'}

""" Dictionary of GCAM basin names (keys) and their abbreviations (values). """
gcam_basin_names_and_abbrevations = {
'Africa_East_Central_Coast': 'AfrCstE', 'Africa_Red_Sea_Gulf_of_Aden_Coast': 'AfrCstNE', 'Africa_North_Interior': 'AfrIntN', 'Congo': 'CongoR',
'Lake_Chad': 'LChad', 'Madagascar': 'Madagascar', 'Nile': 'NileR', 'Rift_Valley': 'RiftValley', 'Shebelli_Juba': 'ShebJubR', 
'Africa_North_West_Coast': 'AfrCstNW', 'Dead_Sea': 'DeadSea', 'Mediterranean_Sea_East_Coast': 'MeditE', 'Mediterranean_South_Coast': 'MeditS', 
'Niger': 'NigerR', 'Sinai_Peninsula': 'SinaiP', 'Africa_Indian_Ocean_Coast': 'AfrCstSE', 'Namibia_Coast': 'AfrCstSW', 
'South_Africa_South_Coast': 'AfrCstS', 'Africa_South_Interior': 'AfrIntS', 'Angola_Coast': 'AngolaCst', 'Gulf_of_Guinea': 'GuineaGulf', 
'Limpopo': 'LimpopoR', 'Orange': 'OrangeR', 'Zambezi': 'ZambeziR', 'Africa_West_Coast': 'AfrCstW', 'Senegal': 'SenegalR', 'Volta': 'VoltaR', 
'South_America_Colorado': 'ArgColoR', 'North_Argentina_South_Atlantic_Coast': 'ArgCstN', 'South_Argentina_South_Atlantic_Coast': 'ArgCstS', 
'North_Chile_Pacific_Coast': 'ChileCstN', 'South_Chile_Pacific_Coast': 'ChileCstS', 'La_Puna_Region': 'LaPuna', 'Mar_Chiquita': 'MarChiq', 
'Negro': 'NegroR', 'Pampas_Region': 'Pampas', 'Central_Patagonia_Highlands': 'Patagonia', 'La_Plata': 'RioLaPlata', 'Salinas_Grandes': 'Salinas', 
'Australia_East_Coast': 'AusCstE', 'Australia_North_Coast': 'AusCstN', 'Australia_South_Coast': 'AusCstS', 'Australia_West_Coast': 'AusCstW', 
'Australia_Interior': 'AusInt', 'Murray_Darling': 'MurrayDrlg', 'New_Zealand': 'NewZealand', 'Tasmania': 'Tasmania', 'Amazon': 'AmazonR', 
'East_Brazil_South_Atlantic_Coast': 'BrzCstE', 'North_Brazil_South_Atlantic_Coast': 'BrzCstN', 'Uruguay_Brazil_South_Atlantic_Coast': 'BrzCstS', 
'Orinoco': 'OrinocoR', 'Parnaiba': 'ParnaibaR', 'Northeast_South_America_South_Atlantic_Coast': 'SAmerCstNE', 'Sao_Francisco': 'SaoFrancR', 
'Tocantins': 'TocantinsR', 'Atlantic_Ocean_Seaboard': 'CanAtl', 'Churchill': 'ChurchillR', 'Fraser': 'FraserR', 'Great_Lakes_Basin': 'GreatLakes', 
'Hudson_Bay_Coast': 'HudsonBay', 'Mackenzie': 'Mackenzie', 'Missouri_River_Basin': 'MissouriR', 'Northwest_Territories': 'NWTerr', 
'Saskatchewan_Nelson': 'NelsonR', 'Pacific_and_Arctic_Coast': 'PacArctic', 'St_Lawrence': 'StLwrncR', 'New_England_Basin': 'UsaCstNE', 
'Pacific_Northwest_Basin': 'UsaPacNW', 'Caribbean': 'Caribbean', 'Southern_Central_America': 'CntAmer', 
'Colombia_Ecuador_Pacific_Coast': 'ColEcuaCst', 'Grijalva_Usumacinta': 'GrijUsuR', 'Caribbean_Coast': 'SAmerCstN', 'Yucatan_Peninsula': 'YucatanP', 
'Amu_Darya': 'AmuDaryaR', 'Amur': 'AmurR', 'Black_Sea_North_Coast': 'BlackSeaN', 'Black_Sea_South_Coast': 'BlackSeaS', 
'Caspian_Sea_East_Coast': 'CaspianE', 'Caspian_Sea_Coast': 'CaspianNE', 'Caspian_Sea_South_West_Coast': 'CaspianSW', 'Gobi_Interior': 'Gobi', 
'Lake_Balkash': 'LBalkash', 'Ob': 'ObR', 'Syr_Darya': 'SyrDaryaR', 'Tarim_Interior': 'Tarim', 'Ural': 'UralR', 'Volga': 'VolgaR', 
'Yenisei': 'YeniseiR', 'Bo_Hai_Korean_Bay_North_Coast': 'BoHai', 'China_Coast': 'ChinaCst', 'Ganges_Bramaputra': 'GangesR', 'Hainan': 'Hainan',
'Hong_(Red_River)': 'Hong', 'Huang_He': 'HuangHeR', 'Indus': 'IndusR', 'Irrawaddy': 'IrrawaddyR', 'Mekong': 'Mekong', 
'Russia_South_East_Coast': 'RusCstSE', 'South_China_Sea_Coast': 'SChinaSea', 'Salween': 'Salween', 'Plateau_of_Tibet_Interior': 'Tibet', 
'Xun_Jiang': 'XunJiang', 'Yangtze': 'Yangtze', 'Ziya_He_Interior': 'ZiyaHe', 'Magdalena': 'MagdalenaR', 
'Adriatic_Sea_Greece_Black_Sea_Coast': 'AdrBlkSea', 'Baltic_Sea_Coast': 'BalticSea', 'Danube': 'DanubeR', 'Daugava': 'DaugavaR', 
'Dniester': 'DniesterR', 'Denmark_Germany_Coast': 'DnkGrmCst', 'Elbe': 'ElbeR', 'Narva': 'NarvaR', 'Neman': 'NemanR', 'Oder': 'OderR', 
'Poland_Coast': 'PolandCst', 'Wisla': 'WislaR', 'Arctic_Ocean_Islands': 'ArcticIsl', 'Douro': 'DouroR', 'Ebro': 'EbroR', 'Ems_Weser': 'EmsWeserR',
'England_and_Wales': 'EngWales', 'Finland': 'Finland', 'France_South_Coast': 'FranceCstS', 'France_West_Coast': 'FranceCstW', 'Gironde': 'Gironde',
'Guadalquivir': 'GuadalqR', 'Guadiana': 'GuadianaR', 'Spain_Portugal_Atlantic_Coast': 'IberiaCst', 'Ireland': 'Ireland', 
'Italy_East_Coast': 'ItalyCstE', 'Italy_West_Coast': 'ItalyCstW', 'Loire': 'LoireR', 'Mediterranean_Sea_Islands': 'MeditIsl', 'Neva': 'NevaR', 
'Po': 'PoR', 'Rhine': 'RhineR', 'Rhone': 'RhoneR', 'Scheldt': 'ScheldtR', 'Scandinavia_North_Coast': 'ScndnvN', 'Scotland': 'Scotland', 
'Seine': 'SeineR', 'Spain_South_and_East_Coast': 'SpainCstSE', 'Sweden': 'Sweden', 'Tagus': 'TagusR', 'Tiber': 'TiberR', 'Dnieper': 'DnieperR', 
'Don': 'DonR', 'Eastern_Jordan_Syria': 'EJrdnSyr', 'Tigris_Euphrates': 'TigrEuphR', 'Iceland': 'Iceland', 'Andaman_Nicobar_Islands': 'AdnNicIsl', 
'Bay_of_Bengal_North_East_Coast': 'BengalBay', 'Yasai': 'BengalW', 'Brahmani': 'BrahmaniR', 'Cauvery': 'CauveryR', 'Godavari': 'GodavariR', 
'India_East_Coast': 'IndCstE', 'India_North_East_Coast': 'IndCstNE', 'India_South_Coast': 'IndCstS', 'India_West_Coast': 'IndCstW', 
'Krishna': 'KrishnaR', 'Mahanadi': 'MahanadiR', 'Mahi': 'MahiR', 'Narmada': 'NarmadaR', 'Pennar': 'PennarR', 'Sabarmati': 'SabarmatiR', 
'Tapti': 'TaptiR', 'North_Borneo_Coast': 'BorneoCstN', 'Fly': 'FlyR', 'Palau_and_East_Indonesia': 'IdnE', 'Irian_Jaya_Coast': 'IrianJaya', 
'Java_Timor': 'JavaTimor', 'Kalimantan': 'Kalimantan', 'Sepik': 'SepikR', 'Sulawesi': 'Sulawesi', 'Sumatra': 'Sumatra', 'Japan': 'Japan', 
'Taiwan': 'Taiwan', 'California_River_Basin': 'California', 'Baja_California': 'MexBaja', 'Mexico_Northwest_Coast': 'MexCstNW', 
'Pacific_Central_Coast': 'MexCstW', 'North_Gulf': 'MexGulf', 'Mexico_Interior': 'MexInt', 'Papaloapan': 'Papaloapan', 'Rio_Balsas': 'RioBalsas', 
'Rio_Grande_River_Basin': 'RioGrande', 'Rio_Lerma': 'RioLerma', 'Rio_Verde': 'RioVerde', 'Isthmus_of_Tehuantepec': 'Tehuantpc', 
'Lower_Colorado_River_Basin': 'UsaColoRS', 'Arabian_Peninsula': 'ArabianP', 'Arabian_Sea_Coast': 'ArabianSea', 'Farahrud': 'FarahrudR', 
'HamuniMashkel': 'HamuMashR', 'Helmand': 'Helmand', 'Central_Iran': 'Iran', 'Persian_Gulf_Coast': 'PersianGulf', 'Red_Sea_East_Coast': 'RedSeaE', 
'Russia_Barents_Sea_Coast': 'BarentsSea', 'Northern_Dvina': 'DvinaRN', 'Kara_Sea_Coast': 'KaraSea', 'Lena': 'LenaR', 
'Siberia_North_Coast': 'SiberiaN', 'Siberia_West_Coast': 'SiberiaW', 'South_Africa_West_Coast': 'AfrCstSSW', 'Peru_Pacific_Coast': 'PeruCst', 
'Sri_Lanka': 'SriLanka', 'North_and_South_Korea': 'Korea', 'Chao_Phraya': 'ChaoPhrR', 'Peninsula_Malaysia': 'MalaysiaP', 
'South_Pacific_Islands': 'NewCaledn', 'Papua_New_Guinea_Coast': 'PapuaCst', 'Philippines': 'Phlppns', 'Sittaung': 'SittaungR', ''
'Solomon_Islands': 'SolomonIsl', 'Gulf_of_Thailand_Coast': 'ThaiGulf', 'Viet_Nam_Coast': 'VietnamCst', 'Arkansas_White_Red_Basin': 'ArkWhtRedR', 
'Great_Basin': 'GreatBasin', 'Hawaii': 'Hawaii', 'Upper_Mississippi_Basin': 'MissppRN', 'Lower_Mississippi_River_Basin': 'MissppRS', 
'Ohio_River_Basin': 'OhioR', 'Tennessee_River_Basin': 'TennR', 'Texas_Gulf_Coast_Basin': 'TexasCst', 'Upper_Colorado_River_Basin': 'UsaColoRN', 
'Mid_Atlantic_Basin': 'UsaCstE', 'South_Atlantic_Gulf_Basin': 'UsaCstSE'
}


def produce_dataframe_for_landtype_group(df, category, category_label, value_label, 
                landtype_groups, mean_or_sum_if_more_than_one_row_in_same_landtype_group, key_columns):
    """ 
    Aggregates the rows of a given Pandas DataFrame that match the specified landtype_group (e.g., crop, forest, pasture, shrub, grass).
    Performs one of four user-specified operations on the rows in the group: mean, sum, area-weighted mean, or area-weighted sum.

    Parameters:
        df: DataFrame containing the data of interest.
        category: String specifying the name of the category of interest.
        category_label: String specifying the label for the appropriate category (e.g., 'sector' or 'landtype').
        value_label: String specifying the label for the column containing the value of interest.
        landtype_groups: Dictionary where the keys are landtype group names and the values are all the landtypes that belong to each group.
        mean_or_sum_if_more_than_one_row_in_same_landtype_group: String that indicates the operation that should be performed on each group.
        key_columns: Columns on which the aggregation (group-by) operation should be performed.

    Returns:
        DataFrame with aggregated rows and the value_label column modified to reflect a mean, sum, area-weighted mean, or area-weighted sum.
    """
    # Get all landtypes for the group and parse the DataFrame to have only the rows that correspond to one of the landtypes in the group.
    landtypes = landtype_groups[category]
    df = df[df[category_label].isin(landtypes)]

    if mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'mean':
        landtypes_in_df = [x for x in landtypes if x in df[category_label].unique()]
        num_landtypes_in_df = len(landtypes_in_df)
        df = df.groupby(key_columns).sum()
        df.loc[:, df.columns != category_label] /= num_landtypes_in_df
        df = df.reset_index()
    elif mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'sum':
        df = df.groupby(key_columns).sum().reset_index()
    elif mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'area_weighted_mean':
        df.loc[:, value_label] = df['area']*df[value_label]
        df = df.groupby(key_columns).sum()
        total_area = df['area']
        df.loc[:, df.columns != category_label] = df.loc[:, df.columns != category_label].div(total_area, axis=0)
        df = df.reset_index()
    elif mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'area_weighted_sum':
        df[value_label] = df['area']*df[value_label]
        df = df.groupby(key_columns).sum().reset_index()
    df[category_label] = category
    return df

def standardize_crop_names(df, columns, mean_or_sum_if_more_than_one_row_for_crop_name='mean'):
    """
    Applies the mappings in gcam_crop_names to produce a common set of crop names in the given Pandas DataFrame.
    An aggregation followed by a mean or sum is performed if there happens to be more than one row that matches a value for the given set of columns.

    Parameters:
        df: DataFrame to modify.
        columns: Columns over which to perform the aggregation (group-by) operation.
        mean_or_sum_if_more_than_one_row_for_crop_name: Specifies whether to calculate a mean or a sum after performing the aggregation.

    Returns:
        DataFrame with the crop names modified so that they belong to the standard common set.
    """
    df = df.replace(gcam_crop_names)
    if mean_or_sum_if_more_than_one_row_for_crop_name == 'mean':
        return df.groupby(columns).mean().reset_index()
    elif mean_or_sum_if_more_than_one_row_for_crop_name == 'sum':
        return df.groupby(columns).sum().reset_index()