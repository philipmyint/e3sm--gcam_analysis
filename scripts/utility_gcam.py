import numpy as np
import pandas as pd

""" GCAM land type groups and aggregations. """
forest_names = ['Forest', 'ProtectedUnmanagedForest', 'UnmanagedForest']
pasture_names = ['Pasture', 'UnmanagedPasture', 'ProtectedUnmanagedPasture']
grassland_names = ['Grassland', 'ProtectedGrassland']
crop_names_original = ['biomass', 'biomassGrass', 'biomassTree', 'Corn', 'CornC4', 'FiberCrop', 'FodderGrass', 'FodderHerb', 'FodderHerbC4', 'Fruits', 
	'FruitsTree', 'Legumes', 'MiscCrop', 'MiscCropC4', 'MiscCropTree', 'NutsSeeds', 'NutsSeedsTree', 'OilCrop', 'OilCropTree', 'OilPalmTree', 
    'OtherGrain', 'OtherGrainC4', 'Rice', 'RootTuber', 'Soybean', 'SugarCrop', 'SugarCropC4', 'Vegetables', 'Wheat'] # Original crop names.
crop_names = ['BioenergyCrop', 'Corn', 'FiberCrop', 'FodderGrass', 'FodderHerb', 'Fruits', 'Legumes', 'MiscCrop', 'NutsSeeds',  
    'OilCrop', 'OilPalm', 'OtherGrain', 'Rice', 'RootTuber', 'Soybean', 'SugarCrop', 'Vegetables', 'Wheat']   # Modified set of crop names.
other_arable_names = ['OtherArableLand']		# This is considered a managed cropland type.
shrubland_names = ['Shrubland', 'ProtectedShrubland']
urban_names = ['UrbanLand']                     # These are constant.
other_names = ['RockIceDesert', 'Tundra']       # These are constant.

""" Unique sector values found in the ag_commodity_prices GCAM output. Used to detect ag-sector DataFrames. """
ag_sectors = ['biomass', 'Corn', 'FiberCrop', 'FodderGrass', 'FodderHerb', 'Forest', 'Fruits', 'Legumes', 'MiscCrop',
              'NutsSeeds', 'OilCrop', 'OilPalm', 'OtherGrain', 'Pasture', 'Rice', 'RootTuber', 'Soybean',
              'SugarCrop', 'UnmanagedLand', 'Vegetables', 'Wheat']

""" Dictionary of GCAM crop name mappings (keys = old names, values = ag_sectors) for adding area to certain dataframes. """
ag_sector_mappings = {'biomass': 'biomass', 'biomassGrass': 'biomass', 'biomassTree': 'biomass', 'CornC4': 'Corn', 
        'FodderHerbC4': 'FodderHerb', 'FruitsTree': 'Fruits', 'MiscCropC4': 'MiscCrop', 'MiscCropTree': 'MiscCrop', 'NutsSeedsTree': 'NutsSeeds',
        'OilCropTree': 'OilCrop', 'OilPalmTree': 'OilPalm', 'OtherGrainC4': 'OtherGrain', 'SugarCropC4': 'SugarCrop',
        'ProtectedUnmanagedForest': 'Forest', 'UnmanagedForest': 'Forest', 'ProtectedUnmanagedPasture': 'Pasture', 'UnmanagedPasture': 'Pasture',}

""" Dictionary of GCAM land type groups and aggregations with modified/standardized crop names. """
gcam_landtype_groups = {'forest': forest_names, 'pasture': pasture_names, 'grass': grassland_names, 'crop': crop_names,
        'other_arable': other_arable_names, 'shrub': shrubland_names, 'urban': urban_names, 'other': other_names}

""" Dictionary of GCAM land type groups and aggregations with original crop names. """
gcam_landtype_groups_original = gcam_landtype_groups.copy()
gcam_landtype_groups_original['crop'] = crop_names_original

""" Dictionary of GCAM crop name mappings (keys = old names, values = new names) to create a standardized set of crop names across all files. """
gcam_crop_mappings = {'biomass': 'BioenergyCrop', 'biomassGrass': 'BioenergyCrop', 'biomassTree': 'BioenergyCrop', 'CornC4': 'Corn', 
        'FodderHerbC4': 'FodderHerb', 'FruitsTree': 'Fruits', 'MiscCropC4': 'MiscCrop', 'MiscCropTree': 'MiscCrop', 'NutsSeedsTree': 'NutsSeeds',
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
'South_Pacific_Islands': 'NewCaledn', 'Papua_New_Guinea_Coast': 'PapuaCst', 'Philippines': 'Phlppns', 'Sittaung': 'SittaungR', 
'Solomon_Islands': 'SolomonIsl', 'Gulf_of_Thailand_Coast': 'ThaiGulf', 'Viet_Nam_Coast': 'VietnamCst', 'Arkansas_White_Red_Basin': 'ArkWhtRedR', 
'Great_Basin': 'GreatBasin', 'Hawaii': 'Hawaii', 'Upper_Mississippi_Basin': 'MissppRN', 'Lower_Mississippi_River_Basin': 'MissppRS', 
'Ohio_River_Basin': 'OhioR', 'Tennessee_River_Basin': 'TennR', 'Texas_Gulf_Coast_Basin': 'TexasCst', 'Upper_Colorado_River_Basin': 'UsaColoRN', 
'Mid_Atlantic_Basin': 'UsaCstE', 'South_Atlantic_Gulf_Basin': 'UsaCstSE'
}

def modify_crop_names(df, columns, mean_or_sum_if_more_than_one_row_for_crop_name='area_weighted_mean'):
    """
    Applies the mappings in gcam_crop_mappings dictionary to produce a modified set of crop names in the given Pandas DataFrame.
    An aggregation followed by a mean or sum is performed if there happens to be more than one row that matches a value for the given set of columns.

    Parameters:
        df: DataFrame to modify.
        columns: Columns over which to perform the aggregation (group-by) operation.
        mean_or_sum_if_more_than_one_row_for_crop_name: Specifies whether to calculate an area-weighted mean, a mean, or a sum after performing the aggregation.

    Returns:
        DataFrame with the crop names modified so that they belong to the modified common set.
    """
    df = df.replace(gcam_crop_mappings)
    # Identify string columns that are not groupby keys; these cannot be aggregated numerically
    # and must be saved and rejoined after groupby operations.
    str_non_key = [c for c in df.columns if c not in columns and df[c].dtype == object]
    str_first = df.groupby(columns)[str_non_key].first().reset_index() if str_non_key else None
    df_numeric = df.drop(columns=str_non_key)
    if mean_or_sum_if_more_than_one_row_for_crop_name == 'mean':
        mean_df = df_numeric.groupby(columns).mean().reset_index()
        # Areas should be summed, not averaged, even in the 'mean' case.
        if 'area' in df.columns:
            area_sum = df_numeric.groupby(columns)['area'].sum().reset_index()
            # Replace/merge the area column in mean_df.
            mean_df = mean_df.drop(columns='area', errors='ignore').merge(area_sum, on=columns)
        if str_first is not None:
            mean_df = mean_df.merge(str_first, on=columns, how='left')
        return mean_df
    elif mean_or_sum_if_more_than_one_row_for_crop_name == 'sum':
        result = df_numeric.groupby(columns).sum().reset_index()
        if str_first is not None:
            result = result.merge(str_first, on=columns, how='left')
        return result
    elif mean_or_sum_if_more_than_one_row_for_crop_name == 'area_weighted_mean':
        # Identify columns to be averaged (these are numeric, non-area columns that are not part of the aggregation operation).
        numeric_cols = df_numeric.select_dtypes('number').columns.tolist()
        numeric_cols = [col for col in numeric_cols if col not in columns + ['area']]

        # Define a function to compute area-weighted mean for each group.
        def weighted_mean(g):
            result = {}
            total_area = g['area'].sum()
            for col in numeric_cols:
                if total_area != 0:
                    result[col] = (g[col] * g['area']).sum() / total_area
                else:
                    result[col] = 0
                    # result[col] = float('nan')  # Avoid division by zero.
            result['area'] = total_area
            return pd.Series(result)

        # Apply the weighted mean function to each group.
        result = df_numeric.groupby(columns).apply(weighted_mean, include_groups=False).reset_index()
        if str_first is not None:
            result = result.merge(str_first, on=columns, how='left')
        return result
    
def produce_dataframe_for_landtype_group(df, category, category_label, value_label, 
                landtype_groups, mean_or_sum_if_more_than_one_row_in_same_landtype_group, key_columns):
    """ 
    Aggregates the rows of a given Pandas DataFrame that match the specified landtype_group (e.g., crop, forest, pasture, shrub, grass).
    Performs one of four user-specified operations on the rows in the group: mean, sum, area-weighted mean, or area-weighted sum.

    Parameters:
        df: DataFrame containing the data of interest.
        category: String specifying the name of the landtype category of interest (e.g., 'forest', 'shrub', 'pasture').
        category_label: String specifying the label for the landtype column in the DataFrame (most likely this will just be 'landtype').
        value_label: String specifying the label for the column containing the value of interest.
        landtype_groups: Dictionary where the keys are landtype group names and the values are all the landtypes that belong to each group.
        mean_or_sum_if_more_than_one_row_in_same_landtype_group: String that indicates the operation that should be performed on each group.
        key_columns: Columns on which the aggregation (group-by) operation should be performed.

    Returns:
        DataFrame with aggregated rows and the value_label column modified to reflect a mean, sum, area-weighted mean, or area-weighted sum.
    """
    # Get all landtypes for the group and filter the DataFrame to keep only the rows that correspond to one of the landtypes in the group.
    landtypes = landtype_groups[category]
    df = df[df[category_label].isin(landtypes)]

    if mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'mean':
        landtypes_in_df = [x for x in landtypes if x in df[category_label].unique()]
        num_landtypes_in_df = len(landtypes_in_df)
        df = df.groupby(key_columns).sum()
        df.loc[:, value_label] /= num_landtypes_in_df
        df = df.reset_index()
    elif mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'sum':
        df = df.groupby(key_columns).sum().reset_index()
    elif mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'area_weighted_mean':
        df.loc[:, value_label] = df['area']*df[value_label]
        df = df.groupby(key_columns).sum()
        total_area = df['area']
        df.loc[:, value_label] = df.loc[:, value_label].div(total_area, axis=0)
        df = df.reset_index()
    elif mean_or_sum_if_more_than_one_row_in_same_landtype_group == 'area_weighted_sum':
        df.loc[:, value_label] = df['area']*df[value_label]
        df = df.groupby(key_columns).sum().reset_index()
    df[category_label] = category
    return df

def add_areas_to_dataframe(df, df_land):
    """
    Adds areas from a processed detailed land allocation DataFrame as an extra 'area' column to df.
    Merges on all columns common to df and df_land, excluding non-key columns ('value', 'area').
    'sector' in df is treated as a placeholder and excluded from merge keys even if present in both.
    If df has a 'market' column and no merge matches are found, falls back to matching region names
    within the market text string.

    Parameters:
        df: DataFrame to add areas to (loaded from a CSV produced by
            gcam_extract_csv_from_xml_or_project_files.R).
        df_land: Processed detailed land allocation DataFrame. Must contain at minimum the columns
                 'scenario', 'region', 'year', and 'value'. Additional key columns (e.g. 'basin',
                 'landtype', and any management columns) are determined by the land allocation file
                 and are used automatically as merge keys where present in both DataFrames.

    Returns:
        DataFrame with added 'area' and 'area_units' columns.
    """
    # Validate that df_land has the minimum columns needed for any merge path.
    required_land_cols = ['scenario', 'region', 'year', 'value']
    missing = [c for c in required_land_cols if c not in df_land.columns]
    if missing:
        print(f"Error: df_land is missing required columns: {missing}. Available: {list(df_land.columns)}")
        return df

    # Rename 'value' -> 'area' and 'units' -> 'area_units' in the land allocation DataFrame so
    # that df's own 'units' column (for the value data) is not used as a merge key and is not
    # overwritten. The land area units are carried through as a separate 'area_units' column.
    rename_map = {}
    if 'value' in df_land.columns:
        rename_map['value'] = 'area'
    if 'units' in df_land.columns:
        rename_map['units'] = 'area_units'
    area_df = df_land.rename(columns=rename_map) if rename_map else df_land.copy()

    # Drop any existing 'area' or 'area_units' columns to avoid duplicates after merge.
    for col in ['area', 'area_units']:
        if col in df.columns:
            df = df.drop(col, axis=1)

    # Determine merge columns: all columns common to df and area_df, excluding non-key columns.
    # 'units' is excluded so df's own units are preserved; area units come through as 'area_units'.
    exclude_from_merge = {'area', 'area_units', 'value', 'units'}
    merge_cols = [c for c in df.columns if c in area_df.columns and c not in exclude_from_merge]

    has_area_units = 'area_units' in area_df.columns

    # If the only merge columns are 'scenario' and 'year' and df has a 'market' column, fall back
    # to matching region names within the market text string for a more specific area lookup.
    # Otherwise merge on all common columns, summing areas for any many-to-one matches.
    if set(merge_cols) == {'scenario', 'year'} and 'market' in df.columns:
        regions = area_df['region'].unique()
        agg_dict = {'area': 'sum'}
        if has_area_units:
            agg_dict['area_units'] = 'first'
        area_by_region = area_df.groupby(['scenario', 'region', 'year']).agg(agg_dict).reset_index()

        def get_area_from_market(row):
            matching = [r for r in regions if r in str(row['market'])]
            if not matching:
                return 0
            matched_region = max(matching, key=len)  # prefer the most specific match
            mask = (
                (area_by_region['scenario'] == row['scenario']) &
                (area_by_region['region'] == matched_region) &
                (area_by_region['year'] == row['year'])
            )
            result = area_by_region.loc[mask, 'area']
            return result.sum() if not result.empty else 0

        df['area'] = df.apply(get_area_from_market, axis=1)
        if has_area_units:
            df['area_units'] = area_df['area_units'].iloc[0] if not area_df.empty else ''
    elif 'landtype' not in df.columns and 'sector' in df.columns and df['sector'].isin(ag_sectors).any():
        # Apply crop name mappings first so that original crop names (e.g. 'CornC4') are standardized
        # (e.g. 'Corn') to match the ag_sector values in df before grouping.
        # Then aggregate df_land areas by (scenario, region, ag_sector, year),
        # rename landtype -> sector, then include sector in the merge keys.
        area_df_mapped = area_df.copy()
        area_df_mapped['landtype'] = area_df_mapped['landtype'].replace(ag_sector_mappings)
        agg_dict = {'area': 'sum'}
        if has_area_units:
            agg_dict['area_units'] = 'first'
        area_by_sector = (area_df_mapped
                          .groupby(['scenario', 'region', 'landtype', 'year'], as_index=False)
                          .agg(agg_dict)
                          .rename(columns={'landtype': 'sector'}))
        sector_merge_cols = merge_cols + ['sector']
        df = df.merge(area_by_sector, on=sector_merge_cols, how='left')
    elif merge_cols:
        agg_dict = {'area': 'sum'}
        if has_area_units:
            agg_dict['area_units'] = 'first'
        area_summed = (area_df[merge_cols + list(agg_dict.keys())]
                       .groupby(merge_cols, as_index=False)
                       .agg(agg_dict))
        df = df.merge(area_summed, on=merge_cols, how='left')

        # Fallback: for rows that still have no area match (NaN), retry with a reduced set
        # of merge keys that excludes any merge_col where the unmatched row has an empty
        # string. This handles cases where df has less detail (e.g. water='') than df_land.
        unmatched = df['area'].isna()
        if unmatched.any():
            empty_col_patterns = (df.loc[unmatched, merge_cols]
                                  .apply(lambda row: tuple(c for c in merge_cols if row[c] == ''), axis=1)
                                  .unique())
            for empty_cols in empty_col_patterns:
                if not empty_cols:
                    continue
                fallback_merge_cols = [c for c in merge_cols if c not in empty_cols]
                if not fallback_merge_cols:
                    continue
                fallback_area = (area_df[fallback_merge_cols + list(agg_dict.keys())]
                                 .groupby(fallback_merge_cols, as_index=False)
                                 .agg(agg_dict))
                pattern_mask = unmatched.copy()
                for c in empty_cols:
                    pattern_mask &= (df[c] == '')
                if not pattern_mask.any():
                    continue
                fill_df = df.loc[pattern_mask, fallback_merge_cols].merge(
                    fallback_area, on=fallback_merge_cols, how='left')
                df.loc[pattern_mask, 'area'] = fill_df['area'].values
                if has_area_units:
                    if 'area_units' not in df.columns:
                        df['area_units'] = None
                    df.loc[pattern_mask, 'area_units'] = fill_df['area_units'].values
                unmatched = df['area'].isna()

    df['area'] = df['area'].fillna(0)
    return df