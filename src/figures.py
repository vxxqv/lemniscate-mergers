"""Render manuscript figures from immutable CSV data and exact model formulas."""
from pathlib import Path
import csv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
DATA=ROOT/'computations/main-003'
TEAL='#1F6F5F';GRAY='#6F6F6F';RED='#B5473C'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,
 'axes.titlesize':10,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,
 'axes.spines.right':False,'axes.linewidth':.65,'lines.linewidth':1.5,
 'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
 'pdf.fonttype':42,'svg.fonttype':'none'})

def save(fig,name):
 # The journal text block scales these 7.1-inch figures down to about
 # 6.15 inches. Keep body labels above 8 pt after that final reduction.
 for item in fig.findobj(match=Text):
  if item.get_visible():item.set_fontsize(max(9.6,item.get_fontsize()))
 for ext in ['pdf','svg','png']:
  # Publish each export atomically so a simultaneous manuscript build cannot
  # open a partially written figure.
  target=OUT/f'{name}.{ext}'
  temporary=OUT/f'.{name}.building.{ext}'
  fig.savefig(temporary,format=ext,dpi=230,bbox_inches='tight',pad_inches=.08)
  if ext=='pdf' and b'%%EOF' not in temporary.read_bytes()[-1024:]:
   raise RuntimeError(f'Incomplete PDF export: {name}')
  os.replace(temporary,target)
 plt.close(fig)

def read(name):return list(csv.DictReader(open(DATA/name)))

def landscape():
 fig=plt.figure(figsize=(7.1,3.55))
 ax=fig.add_axes([.015,.18,.435,.70],projection='3d')
 fig.text(.035,.955,'A  Modulus landscape',fontsize=9.5,va='top')
 fig.text(.60,.955,'B  Merger level sets',fontsize=9.5,va='top')
 xx=np.linspace(-1.45,1.45,130);yy=np.linspace(-.9,.9,110)
 x,y=np.meshgrid(xx,yy);z=np.abs((x+1j*y)**2-1)
 ax.computed_zorder=False
 ax.plot_surface(x,y,np.minimum(z,2.6),cmap=LinearSegmentedColormap.from_list('land',['#e5efeb',TEAL]),linewidth=0,alpha=.92,rasterized=True,zorder=1)
 ax.contour(x,y,z,levels=[1],colors=[RED],linewidths=1.1,zorder=3)
 ax.scatter([-1,1],[0,0],[0,0],c='#111111',s=20,depthshade=False,zorder=5)
 ax.scatter([0],[0],[1],c=RED,s=25,depthshade=False,zorder=6)
 ax.set(xlabel='',ylabel='',zlabel='',zlim=(0,2.6),
        xticks=[-1,0,1],yticks=[-.5,0,.5],zticks=[0,1,2])
 ax.view_init(27,-60);ax.tick_params(labelsize=8,pad=0)
 ax.xaxis.labelpad=2;ax.yaxis.labelpad=2
 fig.text(.17,.045,'Re w',ha='center',va='center',fontsize=9)
 fig.text(.425,.185,'Im w',ha='center',va='center',fontsize=9)
 fig.text(.486,.60,r'$|Q(w)|$',rotation=90,ha='center',va='center',fontsize=9)
 ax2=fig.add_axes([.60,.235,.355,.57]); levels=[.8,1,1.2]
 cs=ax2.contour(x,y,z,levels=levels,colors=[GRAY,RED,TEAL],linewidths=[1,1.5,1.5])
 ax2.clabel(cs,fmt='%0.1f',inline=True,fontsize=8)
 ax2.scatter([-1,1],[0,0],c='#111111',s=24,zorder=5)
 ax2.scatter([0],[0],c=RED,s=24,zorder=5)
 ax2.annotate('saddle',xy=(0,0),xytext=(0,.77),ha='center',fontsize=8,
              arrowprops={'arrowstyle':'-','lw':.6},bbox={'facecolor':'white','edgecolor':'none','pad':1})
 ax2.set(aspect='equal',xlabel='Re w',ylabel='Im w',xticks=[-1,0,1],yticks=[-.5,0,.5])
 ax2.yaxis.set_label_position('right');ax2.yaxis.tick_right()
 save(fig,'landscape')

def scaling():
 rows=read('critical_values.csv');families=['aligned4','opposed2','ring4','broken_ring4']
 names=['Aligned (p = 1)','Opposed (p = 2)','Ring (p = 4)','Broken ring']
 styles=[(GRAY,'o','-'),('#111111','s','--'),(TEAL,'D','-'),(RED,'^',':')]
 fig,axes=plt.subplots(1,3,figsize=(7.1,3.4))
 fig.subplots_adjust(left=.09,right=.99,bottom=.29,top=.87,wspace=.45)
 for ax,key,title in zip(axes,['displacement','frozen_value_error','newton_value_error'],
                          ['A  Displacement','B  Frozen-value error','C  Quadratic-value error']):
  for fam,name,(col,mk,ls) in zip(families,names,styles):
   rr=[r for r in rows if r['family']==fam and r['gamma']=='1']
   rr=sorted(rr,key=lambda r:float(r['epsilon']))
   ax.loglog([float(r['epsilon']) for r in rr],[float(r[key]) for r in rr],color=col,marker=mk,markersize=3,ls=ls,label=name)
  ax.set(xlabel=r'Cluster scale $\varepsilon$')
  ax.set_title(title,fontsize=9,pad=9)
  ax.set_xticks([.02,.1,.5],['0.02','0.1','0.5'])
  ax.grid(axis='y',color='#E8E8E5',linewidth=.5);ax.minorticks_off()
 handles,labels=axes[0].get_legend_handles_labels()
 fig.supylabel('Absolute error',x=.006,fontsize=9)
 fig.legend(handles,labels,fontsize=8,loc='lower center',bbox_to_anchor=(.53,.005),ncol=2,frameon=False)
 save(fig,'scaling')

def hierarchy():
 rows=read('hierarchy_bounds.csv');fams=['aligned4','opposed2','ring4','broken_ring4']
 es=['1/2','1/3','1/5','1/8','1/12','1/20','1/32','1/50']
 fig,axes=plt.subplots(2,1,figsize=(7.1,3.8))
 fig.subplots_adjust(left=.18,right=.985,bottom=.15,top=.90,hspace=.68)
 for ax,key,title in zip(axes,['distance_certifies_hierarchy','moment_certifies_hierarchy'],
                        ['A  Distance-only bracket','B  Moment-aware bracket']):
  a=np.array([[int(next(r for r in rows if r['family']==f and r['epsilon_exact']==e)[key]=='True') for e in es] for f in fams])
  ax.imshow(a,vmin=0,vmax=1,cmap=ListedColormap(['#F1F1EE',TEAL]),aspect='auto')
  for i in range(4):
   for j in range(8):ax.text(j,i,'Y' if a[i,j] else '-',ha='center',va='center',color='white' if a[i,j] else GRAY,fontsize=9)
  ax.set(yticks=range(4),yticklabels=['Aligned','Opposed','Ring','Broken ring'],xticks=range(8),xticklabels=es)
  ax.set_title(title,fontsize=9.5,pad=9)
  ax.tick_params(length=0)
 fig.supxlabel(r'Cluster scale $\varepsilon$',y=.02,fontsize=9)
 save(fig,'hierarchy')

def raster():
 rows=read('raster_trials.csv');ratios=[.125,.25,.5,1,2,4,8]
 fig,axes=plt.subplots(1,2,figsize=(7.1,3.5))
 fig.subplots_adjust(left=.085,right=.985,bottom=.29,top=.87,wspace=.26)
 for ax,conn in zip(axes,[4,8]):
  for stat,col,ls,mk in [('correct',TEAL,'-','o'),('wrong',RED,'--','s'),('roots_unresolved',GRAY,':','^')]:
   vals=[]
   for ratio in ratios:
    rr=[r for r in rows if int(r['adjacency'])==conn and float(r['h_over_width'])==ratio]
    vals.append(sum(r['status']==stat for r in rr)/len(rr))
   ax.plot(ratios,vals,color=col,linestyle=ls,marker=mk,markersize=4,label=stat.replace('_',' '))
  ax.axvline(1/np.sqrt(2),color='#111111',ls=':',lw=1)
  ax.set(xscale='log',ylim=(-.035,1.06))
  ax.set_title(f'{"A" if conn==4 else "B"}  {conn}-neighbor graph',fontsize=9.5,pad=9)
  ax.set_xticks([.125,.5,2,8],['1/8','1/2','2','8']);ax.minorticks_off();ax.grid(axis='y',color='#E8E8E5',lw=.5)
 axes[0].set_ylabel('Fraction of trials')
 handles,labels=axes[0].get_legend_handles_labels()
 fig.legend(handles,labels,frameon=False,fontsize=8,ncol=3,loc='lower center',bbox_to_anchor=(.53,.005))
 fig.supxlabel(r'Grid spacing / feature width, $h/W$',y=.135,fontsize=9)
 save(fig,'raster')

def fields():
 fig,axes=plt.subplots(1,3,figsize=(7.1,2.95))
 fig.subplots_adjust(left=.08,right=.985,bottom=.19,top=.77,wspace=.37)
 roots=np.roots([1,0,-3,.4+.2j]);B=3+1j
 for ax,angs,title in zip(axes,[[0,0,0,0],[0,0,np.pi,np.pi],[0,np.pi/2,np.pi,3*np.pi/2]],['A  Aligned: p = 1\nRoot multiplicity 4','B  Opposed: p = 2\nRoot multiplicity 2 each','C  Ring: p = 4\nSimple external roots']):
  rr=B*np.exp(1j*np.array(angs));ax.add_patch(plt.Circle((0,0),abs(B),fill=False,ec='#B6B6B0',lw=.7))
  ax.scatter(rr.real,rr.imag,marker='x',c=TEAL,s=45,zorder=5)
  ax.scatter(.25*roots.real,.25*roots.imag,c='#111111',s=12,zorder=5)
  ax.set(xlim=(-3.6,3.6),ylim=(-3.6,3.6),aspect='equal')
  ax.set_title(title,fontsize=9,pad=10,linespacing=1.45)
  ax.set_xticks([-3,0,3]);ax.set_yticks([-3,0,3])
 fig.supylabel('Im z',x=.006,fontsize=9)
 fig.supxlabel('Re z',y=.035,fontsize=9)
 save(fig,'external_fields')

if __name__=='__main__':
 landscape();scaling();hierarchy();raster();fields()
