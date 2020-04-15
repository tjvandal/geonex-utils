## 
# This script performs basic tests on GeoNEX L1G data
# 1. Counts the number of (year,day) pairs for each tile saved to .csv
# 2. Generates daytime images composites for each band saved to .jpg files
# 3. Computes statistics mean, standard deviation, and percent_finite from 
#    a random sample of files. Saves statistics in a .csv

import os, sys
import time

from pyhdf.SD import SD, SDC
import glob
import numpy as np
import pandas as pd
import scipy.ndimage
import matplotlib.pyplot as plt

import argparse


def read_geonexl1g_file(file, solar=False):
    '''
    Read a GeoNEX L1G file with pyhdf
    '''
    fp = SD(file, SDC.READ)
    bands = range(1,17)
    arrs = []
    for b in bands:
        b_obj = fp.select("BAND%02i" % b)
        arr = b_obj.get()[:].astype(np.float32)
        fill_value = 32768. # fill value in file doesn't match
        attrs = b_obj.attributes()
        arr[arr == fill_value] = np.nan
        arr *= attrs['Scale_Factor']
        arr += attrs['Offset_Constant']
        arrs.append(arr)
        
    if solar:
        sa = fp.select('Solar_Azimuth').get()[:]
        sz = fp.select('Solar_Zenith').get()[:]
        s = np.concatenate([np.expand_dims(sz, 2), np.expand_dims(sa, 2)], 2)
        s = s * 0.01
        return arrs, s
    return arrs
        
    
def get_tile_days(data_path, tile):
    '''
    Given a tile, return all (year, day) pairs found in data paths
    '''
    h = tile[1:3]
    v = tile[4:6]
    tile_path = os.path.join(data_path, tile)
    years = sorted([int(y) for y in os.listdir(tile_path) if y[0] == '2'])
    info = []
    for y in years:
        year_path = os.path.join(tile_path, str(y))
        days = sorted([int(d) for d in os.listdir(year_path)])
        info += [[tile, y, d] for d in days]
    df = pd.DataFrame(info, columns=['tile', 'year', 'day'])
    return df

def get_year_day_counts(data_path, output_path):
    '''
    Count the number of (year,day) pairs for every tile and save to a .csv
    '''
    tiles = [p for p in os.listdir(data_path) if p[0] == 'h']
    tiles = tiles[:4]
    dfs = []
    for tile in tiles:
        dfs += [get_tile_days(data_path, tile)]
        
    dfs = pd.concat(dfs)
    
    counts = pd.pivot_table(dfs, values='day', index='tile', columns='year', aggfunc=pd.Series.nunique)
    counts_file = os.path.join(output_path, 'tile_year_day_counts.csv')
    counts.to_csv(counts_file)
    print(f"Saved Tile Year Day Counts to {counts_file}")
    return dfs

        
def make_composite_image(data_path, output_path, hrange, vrange, year, day, hour, minute=0):
    '''
    Generate images at 1km for each band over a set of tiles
    Saves images to output_path/*.jpg
    '''
    tile_size = 600    
    composite = np.empty((tile_size*len(vrange), tile_size*len(hrange), 16))
    
    for i, h in enumerate(hrange):
        for j, v in enumerate(vrange):
            ix = i*600
            jy = j*600
        
            tile = f'h{h:02g}v{v:02g}'
            path = os.path.join(data_path, tile, str(year), '%03i' % day, f'*_{year}*_{hour:02g}{minute:02g}*_*.hdf')
            file_list = sorted(glob.glob(path))
            if len(file_list) == 0:
                continue
            file = file_list[0]
            file_data = read_geonexl1g_file(file)
            
            for b, arr in enumerate(file_data):
                scale = tile_size / arr.shape[0] 
                arr_resize = scipy.ndimage.zoom(arr, scale)
                composite[jy:jy+tile_size,ix:ix+tile_size,b] = arr_resize
    
    for b in range(16):
        img = composite[:,:,b]
        to_file = os.path.join(output_path, f'composite_{year}_{day:03g}_{hour:02g}_C{b+1:02g}.jpg')
        plt.imshow(img, cmap='jet')
        plt.title(os.path.basename(file))
        plt.colorbar()
        plt.savefig(to_file)
        plt.close()
        print(f"Saved image to: {to_file}")
    
    
def random_statistics(data_path, output_path, hrange, vrange, samples=10):
    '''
    Randomly samples files from the data_path.
    Estimates mean, standard deviation, and percent finite saved to summary_statistics*.csv
    '''
    tiles = [f'h{h:02g}v{v:02g}' for h in hrange for v in vrange]
    days = get_tile_days(data_path, tiles[0])
    N_days = days.shape[0]
    
    collect_stats = []
    for i in range(samples):
        rand_tile = np.random.choice(tiles, 1)[0]
        rand_day_row = np.random.choice(range(N_days), 1)[0]
        row = days.iloc[rand_day_row]
        rand_hour = np.random.choice(range(24), 1)[0]
        rand_path = os.path.join(data_path, rand_tile, str(row.year), '%03i' % row.day, f'*_{row.year}*_{rand_hour:02g}*_*.hdf')
        rand_files = glob.glob(rand_path)
        if len(rand_files) == 0:
            continue
            
        file = np.random.choice(rand_files, 1)[0]
        sample, solar = read_geonexl1g_file(file, solar=True)
        for j, s in enumerate(sample):
            mean = np.nanmean(s)
            std = np.nanstd(s)
            percent_finite = 100.*np.mean(np.isfinite(s))
            avg_solar = np.nanmean(solar)
            collect_stats.append(dict(band=j+1,mean=mean,std=std,percent_finite=percent_finite,avg_solar=avg_solar))
            
    stats = pd.DataFrame(collect_stats)
    summary = pd.pivot_table(stats, index='band', values=['mean', 'std', 'percent_finite'])
    summary_file = os.path.join(output_path, f'summary_statistics_{samples}samples.csv')
    summary.to_csv(summary_file)
    return summary
            

def main(args):
    assert args.mode in ['all', 'counts', 'images', 'statistics']

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
        
    vrange = list(range(0,20))
    if args.sensor == 'HIMAWARI8':
        hrange = list(range(44,60)) + list(range(0, 3))
        daytime_hour = 2
    elif args.sensor == 'GOES16':
        hrange = list(range(7, 27))
        daytime_hour = 16
    elif args.sensor == 'GOES17':
        hrange = list(range(57, 60)) + list(range(0, 17))
        daytime_hour = 19
        
    if args.mode in ['all', 'counts']:
        tile_year_days = get_year_day_counts(args.data_path, args.output_path)
        print(tile_year_days)
        
    if args.mode in ['all', 'images']:
        t0 = time.time()
        hrange = hrange[:3]
        vrange = vrange[:3]
        make_composite_image(args.data_path, args.output_path, hrange, vrange, args.year, args.dayofyear, daytime_hour)
        print(f"Time to make composite image={time.time()-t0} seconds")

    if args.mode in ['all', 'statistics']:
        stats = random_statistics(args.data_path, args.output_path, hrange, vrange)
        print(stats)
        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', type=str, default='/nex/projects/goesscratch/weile/AHI05_java/output/HDF4/', 
                        help="Pass base directory of GeoNEX L1G data product.")
    parser.add_argument('--sensor', type=str, default='HIMAWARI8', help="Select from [HIMAWARI8,GOES16,GOES17]")
    parser.add_argument('--output_path', type=str, default='./data_statistics/HIMAWARI8', help='Directory to save files')
    parser.add_argument('--mode', type=str, default='all', help='Select a model from [all, counts, images, statistics]')
    parser.add_argument('--year', type=int, default=2018, help='[images] Select year')
    parser.add_argument('--dayofyear', type=int, default=180, help='[images] Select day of year (1-366)')
    parser.add_argument('--samples', type=int, default=100, help='[statisics] Select number of files to sample')

    args = parser.parse_args()
    
