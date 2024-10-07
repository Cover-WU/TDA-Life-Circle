import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt 
import geopandas as gpd

from datetime import timedelta, date
from functools import reduce
from numba import jit

# * The same function as in module: track_visualization.py, LINKAGE.
def activities_projection(data):
    '''
    parameter: 
        data: a record dataset from one individual.
    '''
    # sort the activities by time...
    data = data.copy()
    # data.sort_values(by = 't_start', ignore_index=True, inplace=True)
    # construct the GeoDataFrame for projection
    data_geocode = gpd.GeoDataFrame(data, geometry=gpd.points_from_xy(data.longitude, data.latitude), crs=4326)
    # extract x-y coordination from UTM projection
    data_geocode.to_crs('EPSG:32650', inplace=True)
    data_geocode['X'] = data_geocode.geometry.x
    data_geocode['Y'] = data_geocode.geometry.y
    data = pd.DataFrame(data_geocode.drop(columns='geometry'))
    return data

# * The same function as in module: track_visualization.py, LINKAGE.
def activity_to_flow(data):
    '''
    parameter: 
        data: a activity record set from one individual.
    '''
    # flow constructing...
    t_origin = data.loc[data.index[:-1], 't_end'].reset_index(drop=True)
    x_origin = data.loc[data.index[:-1], 'X'].reset_index(drop=True)
    y_origin = data.loc[data.index[:-1], 'Y'].reset_index(drop=True)
    t_destination = data.loc[data.index[1:], 't_start'].reset_index(drop=True)
    x_destination = data.loc[data.index[1:], 'X'].reset_index(drop=True)
    y_destination = data.loc[data.index[1:], 'Y'].reset_index(drop=True)
    
    flow = pd.DataFrame({'t_origin': t_origin, 't_destination': t_destination,
        'x_origin': x_origin, 'y_origin': y_origin, 'x_destination': x_destination, 'y_destination': y_destination})
    
    # calculate travelling speed
    # the travel speed of flows generated must be slower than the real value
    flow['duration'] = flow['t_destination'] - flow['t_origin']
    flow['x_diff'] = flow['x_origin'] - flow['x_destination']
    flow['y_diff'] = flow['y_origin'] - flow['y_destination']
    flow['distance'] = np.linalg.norm(flow[['x_diff', 'y_diff']], axis=1)
    flow['speed'] = 3.6 * flow['distance'] / (flow['duration'] / np.timedelta64(1, 's'))
    return flow


def travel_chain_cut(stay_parti: pd.DataFrame):
    '''
    parameter:
        a slice of stay records from one individual.
    '''

    stay_parti = stay_parti.sort_values(by=['t_start', 't_end'])
    move_parti = activity_to_flow(stay_parti) 
    # seperate the move dataset into two parts
    the_same_day = move_parti.t_origin.apply(lambda x : x.day) == move_parti.t_destination.apply(lambda x : x.day)
    move_within_same_day, move_across_day = move_parti[the_same_day], move_parti[~the_same_day]

    # screen out invalid move...
    travel_track_part1 = move_within_same_day.loc[(move_within_same_day.distance > 0) | (move_within_same_day.duration > pd.Timedelta(15, 'm'))]
    travel_track_part2 = move_across_day.loc[(move_across_day.distance > 400) & (move_across_day.t_destination.dt.date - move_across_day.t_origin.dt.date == timedelta(1))]
    travel_track = pd.concat([travel_track_part1, travel_track_part2])
    # ... and just sort it.
    travel_track.sort_values(by = 't_origin', ignore_index=True, inplace=True)
    travel_track['cutting'] = range(1, 1+len(travel_track))

    # cross the valid move with all stays thus cutting the stay series into multiple parts.
    stay_track = stay_parti.copy()
    travel_track.rename(columns={'t_origin':'t_start', 't_destination':'t_end'}, inplace=True)
    track_demo = pd.concat([stay_track, travel_track], ignore_index=True)
    track_demo.sort_values(by=['t_start', 't_end'], ignore_index=True, inplace=True)
    track_demo['cutting'] = track_demo.cutting.fillna(method='pad').fillna(0)
    # track_demo.info()

    # only keep the stay records
    track_stay = track_demo.loc[~track_demo.longitude.isna()].copy()
    # track_stay.drop(columns=track_stay.columns[8:16], inplace=True)
    track_stay = track_stay.dropna(axis=1, how='all') \
        .reindex(columns=['pid', 'cutting', 't_start', 't_end', 'poi_id', 'longitude', 'latitude', 'X', 'Y', 'ptype'])
    return track_stay


# for each cut (subDataFrame), connect the activities.
@jit(nopython=True)
def stay_connection(arr: np.ndarray):
    arr = arr.copy()
    arr_dt = arr[:, 1:3].copy()
    arr_main = arr[:, 3:].copy()
    if arr.shape[0] <= 1:
        return arr
    else:
        date_slice_dt = np.sort(np.unique(np.floor_divide(arr_dt.flatten(), int(1e9 * 3600 * 24))))
        date_slice_pd = np.arange(date_slice_dt.min(), date_slice_dt.max() + 1)
        if len(date_slice_dt) == len(date_slice_pd):
            arr[0,2] = arr[-1,2]
            return arr[0:1, :]
        else:
            arr_current = arr[0:0, :]
            missing_cut = np.where(np.diff(date_slice_dt) > 1)[0] + 1
            rec_start = np.roll(np.append(missing_cut, 0), 1)
            rec_end = np.append(missing_cut, arr.shape[0])
            # rec_slicing = np.column_stack((rec_start, rec_end))
            for id in range(len(missing_cut) + 1):
                r0, r = rec_start[id], rec_end[id]
                # arr_dt_sect = arr_dt[r0:r, :]; arr_main_sect = arr_main[r0:r, :]
                arr_sect = arr[r0:r, :]
                arr_sect[0,2] = arr_sect[-1,2]
                arr_current = np.vstack((arr_current, arr_sect[0:1, :]))
            return arr_current


def travel_chain_one_commune(stay_one_commune: pd.DataFrame):
    stay_one_commune = activities_projection(stay_one_commune)
    stay_one_commune = stay_one_commune.sort_values(by=['pid', 't_start', 't_end'])
    stay_prepared = stay_one_commune.groupby(['pid'], as_index=False).apply(travel_chain_cut).reset_index(drop=True)
    
    # First run to compile
    stay_prepared['grpid'] = stay_prepared.groupby(['pid', 'cutting']).ngroup()
    grp_keys = stay_prepared[['pid', 'cutting', 'grpid']].drop_duplicates()
    stay_fnb = stay_prepared[['grpid', 't_start', 't_end', 'poi_id', 'longitude', 'latitude', 'X', 'Y', 'ptype']].copy()
    stay_fnb['t_start'] = stay_fnb['t_start'].astype('int64');  stay_fnb['t_end'] = stay_fnb['t_end'].astype('int64')
    resarr = stay_fnb.groupby('grpid', as_index=False).apply(lambda df: stay_connection(df.to_numpy()))
    resarr = np.vstack(resarr)
    stay_res = pd.DataFrame(data=resarr, columns=['grpid', 't_start', 't_end', 'poi_id', 'longitude', 'latitude', 'X', 'Y', 'ptype'])
    stay_res = pd.merge(grp_keys, stay_res, on='grpid').drop(columns='grpid')
    stay_res = stay_res.assign(t_start = pd.to_datetime(stay_res.t_start), t_end = pd.to_datetime(stay_res.t_end))
    return stay_res

def stay_connection_wrapper(df: pd.DataFrame):
    return stay_connection(df.to_numpy())