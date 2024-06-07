import os
import geopandas as gpd
from shapely.geometry import Polygon

class vec:
    def __init__(self, input_vec, keep_geom_only=True, range=None, init_area=True, **kwargs) -> None:
        
        if isinstance(input_vec, gpd.GeoDataFrame):
            self.gdf = input_vec
        elif isinstance(input_vec, str) and os.path.isfile(input_vec):
            self.gdf = gpd.read_file(input_vec)
        else:
            raise NotImplementedError

        self.N = self.gdf.shape[0]
        
        if keep_geom_only:
            self._eliminate_field()
        
        self._range(range, **kwargs)
        
        if init_area:
            self.area()
    
    def intersection(self, vec2, keep_geom_type=True):
        return self.gdf.overlay(vec2.gdf, how='intersection', keep_geom_type=keep_geom_type)

    def area(self, record_name='area'):
        self.gdf[record_name] = self.gdf.area
    
    def _eliminate_field(self):
        self.gdf = gpd.GeoDataFrame(
            {'gid': [e for e in range(1, self.N + 1)], 'geometry': self.gdf['geometry']}, crs=self.gdf.crs
        )

    def _range(self, range, **kwargs):
        if range is None:
            # self.range = self.gdf.total_bounds
            self.range = Polygon.from_bounds(*self.gdf.total_bounds)
        elif isinstance(range, list):
            raise NotImplementedError
        elif isinstance(range, str):
            _range_gdf = gpd.read_file(range)
            if 'range_trans' in kwargs.keys():
                self.range = Polygon.from_bounds(*_range_gdf.total_bounds)
            else:
                self.range = _range_gdf['geometry'][0]
                
        elif isinstance(range, Polygon):
            self.range = range
            
            # a = gpd.GeoDataFrame({
            #     'geometry': [self.range]
            # }, crs = self.gdf.crs)
            
            # a.to_file('./bound.shp')
            
    
    def __repr__(self) -> str:
        return f'Vector data with {self.N} elements.'
        

class accuracy:
    def __init__(self, gt_path, pred_path):
        self.gt_path        = gt_path
        self.pred_path      = pred_path
    
    def evaluate(self, metrics=[], type=None):
        pass
    
    def accumulate(self, others, metrics=[], type=None):
        pass
    
class vec2vec_acc(accuracy):
    '''
    This class conducts the object-level assessment
    '''
    AVAILABLE_METRIC = ['precision', 'recall', 'accuracy', 'F1']
    def __init__(self, gt_path:str, pred_path:str):
        super().__init__(gt_path, pred_path)
        
        self.pred    = vec(self.pred_path)
        self.gt      = vec(self.gt_path)
        
        self.seg_metric = dict()
        self.det_metric = dict()
    
    def evaluate(self, metrics=[], type='det', iou_threshold=None):
        '''
        @metircs: see AVAILABLE_METRIC
        @tyep   : det, seg and ouseg for alternating 
            In all the mode, the polygons having the max iou to each gt polygons 
            will be filtered out
            
            In det mode,  the amount of the filtered out polygon is TP. Then:
                recall      = TP / gt.N
                precision   = TP  pred.N
            
            In seg mode, TP is derived by the total area filtered out polygon.
            and the FP and TN are the residuals between TP and total area of 
            pred and gt respectively.
            as for FN, there should be an overall area to connstrain the "False"
            area, so the range will be derived in `class vec`. TODO: if there is 
            no range file or vector claimed in the initialization of the `vec`, 
            the range will be derived by the total boundary of the vec. Then:
                recall      = TP / pred.area
                precesion   = TP /  
            
        '''
        
        if type == 'det':
            metrics = v2v_det(self.pred, self.gt, iou_threshold=iou_threshold)
        elif type == 'seg':
            metrics = v2v_seg(self.pred, self.gt, iou_threshold=iou_threshold)
        elif type == 'ouseg':
            metrics = v2v_ouseg(self.pred, self.gt)
        elif type == 'PoLiS':
            raise NotImplementedError
        
        
        return metrics
                

def v2v_det(pred:vec, gt:vec, iou_threshold=None):
    assert iou_threshold is not None
    intersec = pred.intersection(gt)
    
    intersec['inter_area'] = intersec.area
    intersec['iou'] = intersec['inter_area'] / (intersec['area_1'] + intersec['area_2'] - 2* intersec['inter_area'])
    
    filtered = intersec[intersec['iou'] >= iou_threshold]
    
    filtered = filtered.loc[filtered.groupby('gid_2')['iou'].idxmax()]
    TP = filtered.shape[0]
        
    precision = TP / pred.N
    recall = TP / gt.N
    
    return {'precision': precision, 'recall': recall}

def v2v_seg(pred:vec, gt:vec, iou_threshold=None):
    assert iou_threshold is not None
    intersec = pred.intersection(gt)
    
    intersec['inter_area'] = intersec.area
    intersec['iou'] = intersec['inter_area'] / (intersec['area_1'] + intersec['area_2'] - 2* intersec['inter_area'])
    
    filtered = intersec[intersec['iou'] >= iou_threshold]
    
    filtered = filtered.loc[filtered.groupby('gid_2')['iou'].idxmax()]
    
    TP = filtered['inter_area'].sum()
        
    _Positive = pred.gdf['area'].sum()
    FP = _Positive - TP
    
    _True = gt.gdf['area'].sum()
    FN = _True - TP
    
    TN = pred.range.union(gt.range).area
    TN -= (TN + FP - 2* TP)
    
    precision = float(TP/ (TP + FP))
    recall = float(TP/ (TP+ FN))
    
    return {'precision': precision, 'recall': recall}
        
def v2v_ouseg(pred:vec, gt:vec):
    
    raise NotImplementedError
    
    intersec = pred.intersection(gt)
    
    intersec['inter_area'] = intersec.area
    intersec['iou'] = intersec['inter_area'] / (intersec['area_1'] + intersec['area_2'] - 2* intersec['inter_area'])
    
    filtered = intersec.loc[intersec.groupby('gid_1')['iou'].idxmax()]
    
    print(filtered)
    
    filtered_pred_ids = filtered['gid_1']
    
    return {'GOC': None, 'GUC': None, 'GTC': None}
    


if __name__ == '__main__':
    pred_shp_path = r''
    
    gt_shp_path = r''
    pass

    
    acc = vec2vec_acc(gt_shp_path, pred_shp_path)
    
    print(acc.evaluate(type='ouseg', iou_threshold=0.5))
    
    