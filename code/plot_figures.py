"""Reproduce the revised classification and monthly-flux figures from source data.

No model fitting, data selection, or changes to correction families are performed.
The default command reproduces Figure 5. The optional spatial supplement enables
Figure S1; monthly tower pairs enable Figure S5. Archived temperature and spatial
arrays enable the descriptive Figure 1 summaries. All plotting uses matplotlib.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm
from matplotlib.ticker import MaxNLocator, FuncFormatter, NullFormatter

W = 170 / 25.4
PHASES = ['Early', 'Middle', 'Late']
CODES = [1,2,3,4,5,6,7,8,9,10,11,12,14]
VEGETATION = ['Evergreen needleleaf','Evergreen broadleaf','Deciduous needleleaf',
              'Deciduous broadleaf','Mixed forest','Closed shrubland','Open shrubland',
              'Woody savanna','Savanna','Grassland','Wetland','Cropland','Crop–natural mosaic']
MOISTURE = ['Arid', 'Semi-arid', 'Dry subhumid', 'Humid']
PALETTE = ['#185E45','#32854A','#639A6B','#91BA72','#436D61','#935F40',
           '#C19A72','#8F9953','#C6C574','#D9D7A7','#548C9A','#E0AD4A','#B57765']
MOISTURE_COLORS = ['#B7784C','#DBC18A','#A7BE9C','#548D9A']
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
    'font.size':9.5, 'axes.labelsize':9.5, 'axes.titlesize':10,
    'legend.fontsize':8.8, 'xtick.labelsize':8.8, 'ytick.labelsize':8.8,
    'pdf.fonttype':42, 'svg.fonttype':'none', 'axes.spines.top':False,
    'axes.spines.right':False, 'axes.linewidth':.65, 'lines.linewidth':1.15})


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def label(ax, letter, x=-.10, y=1.035):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10.5, fontweight='bold')


def export(fig, key, output, geometry):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bb = fig.get_tightbbox(renderer)
    w, h = fig.get_size_inches()
    margins = [bb.x0*25.4,(w-bb.x1)*25.4,bb.y0*25.4,(h-bb.y1)*25.4]
    assert min(margins) > 0, (key, margins)
    geometry.append(dict(figure=key, width_mm=w*25.4, height_mm=h*25.4,
                         left_mm=margins[0],right_mm=margins[1],
                         bottom_mm=margins[2],top_mm=margins[3]))
    d = output / key
    d.mkdir(parents=True, exist_ok=True)
    fig.savefig(d/'figure.png', dpi=300, facecolor='white')
    fig.savefig(d/'figure.pdf', dpi=300, facecolor='white')
    fig.savefig(d/'figure.svg', dpi=300, facecolor='white')
    fig.savefig(d/'figure.tiff', dpi=600, facecolor='white',
                pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)
    print(key, flush=True)


def hydro(source, output, geometry):
    all_models = pd.read_csv(source/'new_thermal_analyses.csv')
    estimates = all_models.query("outcome=='WUE' and term=='Thermal_slope'")
    support = pd.read_csv(source/'model_support.csv').fillna('')
    fig, axs = plt.subplots(1,3,figsize=(W,5.25),sharey=True)
    fig.subplots_adjust(left=.245,right=.985,bottom=.31,top=.93,wspace=.085)
    cm = mpl.colormaps['RdBu_r'].copy()
    cm.set_bad('#e6e6e6')
    values = []
    for ax, phase, letter in zip(axs, PHASES, 'abc'):
        mat = np.full((13,4),np.nan)
        for i, code in enumerate(CODES):
            for j in range(4):
                key = f'igbp{code}_moist{j}'
                s = support[(support.analysis==key)&(support.phase==phase)]
                z = estimates[(estimates.analysis==key)&(estimates.phase==phase)]
                assert len(s)==1 and len(z)<=1
                s = s.iloc[0]
                if len(z):
                    z = z.iloc[0]
                    mat[i,j] = z.beta*100
                    if z.q_BH<.05:
                        ax.plot(j,i,'o',color='black',ms=2.7)
                else:
                    assert s.blocks<8 and '8 occupied' in s.reason
                    ax.plot(j,i,'x',color='#777777',ms=4,mew=.8)
                values.append(dict(analysis=key,phase=phase,display_value=mat[i,j],
                                   moisture=MOISTURE[j],n_clusters=s.blocks))
        im = ax.imshow(mat,vmin=-5,vmax=5,cmap=cm,aspect='auto',zorder=-1)
        ax.set_title(phase,pad=8)
        ax.set_xticks(range(4),['Arid','Semi-arid','Dry subhumid','Humid'],
                      rotation=55,ha='right',rotation_mode='anchor')
        ax.set_yticks(range(13),VEGETATION)
        ax.tick_params(length=0,pad=4)
        label(ax,letter)
    cax = fig.add_axes([.30,.09,.65,.025])
    fig.colorbar(im,cax=cax,orientation='horizontal',extend='both',ticks=[-5,0,5],
                 label='WUE thermal slope (100β, °C⁻¹)')
    pd.DataFrame(values).to_csv(output/'Figure_5_plot_values.csv',index=False)
    export(fig,'Figure_5',output,geometry)


def stratification(source, gis_path, output, geometry):
    z = np.load(gis_path, allow_pickle=False)
    s = pd.read_csv(source/'GIS_stratification_support.csv')
    fig = plt.figure(figsize=(W,7.0))
    a = fig.add_axes([.09,.77,.56,.17])
    b = fig.add_axes([.09,.49,.56,.17])
    c = fig.add_axes([.335,.09,.475,.29])
    for ax,field,classes,colors,title in [
        (a,'igbp',CODES,PALETTE,'Vegetation classes'),
        (b,'aridity',list(range(4)),MOISTURE_COLORS,'Moisture backgrounds')]:
        arr = np.full(1400*7200,np.nan)
        v = np.full(len(z[field]),len(classes),float)
        for k, code in enumerate(classes):
            v[z[field]==code] = k
        arr[z['flat_id']] = v
        cm = ListedColormap(colors+['#A7A7A7'])
        cm.set_bad('white')
        ax.imshow(arr.reshape(1400,7200),extent=[-180,180,20,90],origin='upper',
                  aspect='auto',cmap=cm,interpolation='nearest',
                  norm=BoundaryNorm(np.arange(len(classes)+2)-.5,len(classes)+1),
                  rasterized=True)
        ax.set(xlim=(-180,180),ylim=(20,90),title=title,ylabel='Latitude (°N)')
        ax.set_xticks([-180,-90,0,90,180])
        ax.set_xlabel('Longitude (°)')
    a.legend([Patch(color=c) for c in PALETTE]+[Patch(color='#A7A7A7')],
             VEGETATION+['Other / unknown'],loc='upper left',bbox_to_anchor=(1.06,1.18),
             frameon=False,fontsize=8.5,labelspacing=.37,handlelength=.9,borderaxespad=0)
    b.legend([Patch(color=c) for c in MOISTURE_COLORS]+[Patch(color='#A7A7A7')],
             MOISTURE+['Unavailable'],loc='upper left',bbox_to_anchor=(1.04,.85),
             frameon=False,fontsize=8.8,labelspacing=.6,handlelength=.9,borderaxespad=0)
    count, clusters = np.zeros((13,4)),np.zeros((13,4))
    for j, veg in enumerate(CODES):
        for k in range(4):
            v = s[(s.igbp==veg)&(s.aridity_group==k)&(s.phase=='Early')]
            assert len(v)==1
            v=v.iloc[0]
            count[j,k],clusters[j,k]=v.cell_years,v.n_clusters
    im=c.imshow(count,cmap='YlGnBu',norm=LogNorm(vmin=count[count>0].min(),vmax=count.max()),aspect='auto')
    for j in range(13):
        for k in range(4):
            c.text(k,j,f'{int(clusters[j,k])}'+(' ×' if clusters[j,k]<8 else ''),
                   ha='center',va='center',fontsize=8.5,
                   color='white' if im.norm(count[j,k])>.6 else '#222222')
    c.set_yticks(range(13),VEGETATION)
    c.set_xticks(range(4),['Arid','Semi-arid','Dry\nsubhumid','Humid'])
    c.tick_params(length=0,pad=4)
    c.set_title('Sample support: spatial clusters',pad=9,fontsize=9.5)
    c.set_xlabel('Moisture background',labelpad=6)
    cax=fig.add_axes([.855,.12,.022,.23])
    fig.colorbar(im,cax=cax)
    cax.yaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:,.0f}'))
    cax.yaxis.set_minor_formatter(NullFormatter())
    cax.set_title('Cell-years\n(log scale)',fontsize=8.5,pad=8)
    for ax, letter in zip((a,b,c),'abc'):
        label(ax,letter)
    export(fig,'Figure_S1',output,geometry)


def monthly(source,pair_path,output,geometry):
    pairs=pd.read_csv(pair_path)
    metrics=pd.read_csv(source/'tower_metrics.csv')
    fig,axs=plt.subplots(1,2,figsize=(W,3.8))
    fig.subplots_adjust(left=.115,right=.98,bottom=.18,top=.90,wspace=.40)
    records=[]
    for ax,endpoint,letter,color in zip(axs,['GPP','ET'],'ab',['#216B45','#255D86']):
        z=pairs[pairs.endpoint==endpoint]
        ref=metrics[(metrics.analysis==endpoint)&(metrics.scale=='Monthly')].iloc[0]
        assert np.isfinite(z[['tower','product']]).all().all()
        x,y=z.tower.to_numpy(),z['product'].to_numpy()
        weights=1/z.groupby('site').tower.transform('count').to_numpy()
        mx,my=np.average(x,weights=weights),np.average(y,weights=weights)
        corr=np.average((x-mx)*(y-my),weights=weights)/np.sqrt(
            np.average((x-mx)**2,weights=weights)*np.average((y-my)**2,weights=weights))
        rmse=np.sqrt(np.average((y-x)**2,weights=weights))
        assert len(z)==ref.pairs and z.site.nunique()==ref.sites
        assert abs(corr-ref.r)<1e-10 and abs(rmse-ref.rmse)<1e-10
        records.append(dict(endpoint=endpoint,site_months=len(z),sites=z.site.nunique(),
                            correlation=corr,rmse=rmse,r_difference=corr-ref.r,
                            rmse_difference=rmse-ref.rmse))
        ax.scatter(x,y,s=3,alpha=.25,color=color,linewidths=0,rasterized=True)
        low,high=min(x.min(),y.min()),max(x.max(),y.max())
        pad=(high-low)*.04
        limits=[low-pad,high+pad]
        ax.plot(limits,limits,'--',color='.45',lw=.8)
        ax.set(xlim=limits,ylim=limits,aspect='equal')
        unit='g C m⁻² d⁻¹' if endpoint=='GPP' else 'mm d⁻¹'
        ax.set_xlabel(f'Tower {endpoint} ({unit})')
        ax.set_ylabel(f'{endpoint} ({unit})')
        ax.set_title(endpoint,pad=8,fontsize=10)
        ax.text(.035,.97,f'r = {ref.r:.3f}\n95% CI {ref.r_lo:.3f}–{ref.r_hi:.3f}\nRMSE = {ref.rmse:.2f}',
                 transform=ax.transAxes,va='top',fontsize=8.8)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        label(ax,letter,y=1.04)
    pd.DataFrame(records).to_csv(output/'monthly_metric_identity.csv',index=False)
    export(fig,'Figure_S5',output,geometry)


def context_summary(source,context_path,metadata_path,gis_path,output):
    context=np.load(context_path)
    meta=np.load(metadata_path)
    retained=np.load(gis_path,allow_pickle=False)['flat_id']
    assert np.array_equal(context['cells'],meta['cells'])
    assert np.allclose(meta['weight'],np.cos(np.deg2rad(meta['lat'])),atol=0,rtol=0)
    rows,profile=[],[]
    for k,phase in enumerate(PHASES):
        values=context['mean'][:,k]
        ok=np.isfinite(values)
        assert np.array_equal(np.sort(context['cells'][ok]),np.sort(retained))
        rows.append(dict(phase=phase,area_weighted_mean_degC=np.average(values[ok],weights=meta['weight'][ok]),
                         area_fraction_soil_warmer_pct=100*np.average(values[ok]>0,weights=meta['weight'][ok]),
                         n_cells=int(ok.sum())))
        for latitude in np.arange(22.5,90,5):
            band=ok&(meta['lat']>=latitude-2.5)&(meta['lat']<latitude+2.5)
            value=np.average(values[band],weights=meta['weight'][band]) if band.any() else np.nan
            profile.append(dict(phase=phase,latitude=latitude,mean=value))
    df=pd.DataFrame(profile)
    ref=pd.read_csv(source/'background_latitude_profile.csv')
    assert np.array_equal(df[['phase','latitude']],ref[['phase','latitude']])
    assert np.allclose(df['mean'],ref['mean'],equal_nan=True,atol=1e-12,rtol=1e-12)
    pd.DataFrame(rows).to_csv(output/'Figure_1_thermal_summary.csv',index=False)
    df.to_csv(output/'Figure_1_latitude_profile.csv',index=False)
    (output/'Figure_1_summary_definition.json').write_text(json.dumps({
        'spatial_summary':'Weighted mean of the existing valid-year cell means; weights are cos(latitude).',
        'cell_mean':'Within-cell mean over eligible years, separately for each phenological phase.',
        'units':'degrees Celsius', 'statistic_type':'Descriptive, not a regression estimate or significance test.',
        'latitude_profile_identity_tolerance':{'atol':1e-12,'rtol':1e-12},
        'retained_footprint_identity':True,
        'source_hashes':{'thermal_context.npz':sha256(context_path),
                         'cell_metadata.npz':sha256(metadata_path),
                         'background_latitude_profile.csv':sha256(source/'background_latitude_profile.csv')}
    },indent=2),encoding='utf8')


def main():
    ap=argparse.ArgumentParser()
    package=Path(__file__).resolve().parents[1]
    ap.add_argument('--source',type=Path,default=package/'data/source_tables')
    ap.add_argument('--monthly-pairs',type=Path)
    ap.add_argument('--gis-classes',type=Path,
                    help='Optional GIS_retained_classes.npz for Figure S1')
    ap.add_argument('--thermal-context',type=Path)
    ap.add_argument('--cell-metadata',type=Path)
    ap.add_argument('--output',type=Path,default=package/'output/figures')
    args=ap.parse_args()
    if (args.thermal_context is None)!=(args.cell_metadata is None):
        ap.error('--thermal-context and --cell-metadata must be provided together')
    if args.thermal_context is not None and args.gis_classes is None:
        ap.error('--thermal-context also requires --gis-classes to check the retained footprint')
    for path in (args.monthly_pairs,args.gis_classes,args.thermal_context,args.cell_metadata):
        if path is not None and not path.is_file():
            ap.error(f'Input file not found: {path}')
    args.output.mkdir(parents=True,exist_ok=True)
    geometry=[]
    if args.thermal_context is not None:
        context_summary(args.source,args.thermal_context,args.cell_metadata,args.gis_classes,args.output)
    hydro(args.source,args.output,geometry)
    if args.gis_classes is not None:
        stratification(args.source,args.gis_classes,args.output,geometry)
    if args.monthly_pairs is not None:
        monthly(args.source,args.monthly_pairs,args.output,geometry)
    pd.DataFrame(geometry).to_csv(args.output/'geometry.csv',index=False)


if __name__=='__main__':
    main()
