% mphopen('E:/brandon_learning/livelink_learning2.mph')
if false
    import com.comsol.model.*
    import com.comsol.model.util.*;
end
if false && ~exist('mymodel','var')
    mymodel = mphload('E:/brandon_learning/livelink_learning.mph','histmodel');
    mphlaunch histmodel
end

%%
%380 microns is width of slide (at 40X)
length_per_pixel=scale40x();



% Create new model and geometry
% histology_model =ModelUtil.create('myModel');
% g1 = histology_model.geom.create('g1',2);

% Convert histology image into a model and geometry
load('histology_image.mat','hist_image')
% hist_geometry = mphimage2geom(hist_image,0.5,'minarea',200,'scale',length_per_pixel,'rectangle','on');
% hist_geometry = mphimage2geom(hist_image,0.5,'minarea',200,'scale',length_per_pixel,'modeltag','alveoli_model','rtol',1e-3);
hist_geometry = mphimage2geom(hist_image,0.5,'minarea',200,'scale',length_per_pixel,'modeltag','alveoli_model','rtol',1e-1);
figure(3); mphgeom(hist_geometry)

% Get tags of existing model and geometry
geometry_tags = hist_geometry.geom.tags();
geom1=hist_geometry.geom(geometry_tags(1));
geom1_features = geom1.feature().tags;

%preallocate before loop
[xmin,ymin,xmax,ymax]=deal(0);

% Get bounds of current domain
nn_ = 1:(length(geom1_features)-1);



for nn = nn_
    % Get the features of the old geometry
    xx = hist_geometry.geom(geometry_tags(1)).feature(geom1_features(nn)).getDoubleArray('x');
    yy = hist_geometry.geom(geometry_tags(1)).feature(geom1_features(nn)).getDoubleArray('y');
    
    xmin = min([xmin; xx(:)]);    
    ymin = min([ymin; yy(:)]);    
    xmax = max([xmax; xx(:)]);    
    ymax = max([ymax; yy(:)]);
        
%     c(nn) = g1.feature.create(sprintf('c%04.0f',nn),'Polygon');
%     c(nn).set('x',xx*length_per_pixel)
%     c(nn).set('y',yy*length_per_pixel)
%     hist_geometry.geom(geometry_tags(1)).feature(geom1_features(nn)).set('x',xx*length_per_pixel);
%     hist_geometry.geom(geometry_tags(1)).feature(geom1_features(nn)).set('y',xx*length_per_pixel);


end

geom1.feature.create('wetspace','Rectangle')
wetspace = hist_geometry.geom(geometry_tags(1)).feature('wetspace');
wetspace.set('y',-ymax)
wetspace.set('ly',2*(ymax-ymin))
wetspace.set('lx',(xmax-xmin))
% START HERE, TRY TO CHANGE COORDINATES BY MULTIPLYING BY LENGTH PER PIXEL
% geom1.feature('curve1')

% Things to run once
if false
    % Add transient pressure
    hist_geometry.physics.create('actd', 'TransientPressureAcoustics', 'geom1');
    hist_geometry.study.create('std1');
    hist_geometry.study('std1').feature.create('time', 'Transient');
    
    
    % Create functions
    hist_geometry.func.create('an1', 'Analytic');
    hist_geometry.func('an1').set('funcname', 'ty_waveform');
    hist_geometry.func('an1').set('expr', 'ty_envelop(y)*sin(2*pi*freq*y)');
    hist_geometry.func('an1').set('args', 'y');

    
    hist_geometry.func.create('rect1', 'Rectangle');
    hist_geometry.func('rect1').set('funcname', 'ty_envelop');
    hist_geometry.func('rect1').set('lower', 't0_ty');
    hist_geometry.func('rect1').set('upper', 't0_ty + dt_ty');
    hist_geometry.func('rect1').set('smooth', '0.1*dt_ty');
    
    hist_geometry.physics('actd').feature.create('pr1', 'Pressure', 1);
    hist_geometry.physics('actd').feature('pr1').set('p0', 1, 'p_atm*ty_envelop');
    hist_geometry.physics('actd').feature('pr1').set('p0', 1, 'p_atm*ty_envelope');
    hist_geometry.physics('actd').feature('pr1').set('p0', 1, 'p_atm');
%     model.physics('actd').feature('pr2').set('p0', 1, 'ty_waveform2(t[Hz])');
%     model.physics('actd').feature('pr2').selection.set([2]);
    
%     model.physics('acpr').feature.create('pr1', 'Pressure', 1);
%     model.physics('acpr').feature('pr1').set('p0', 1, 'p_sound*ty_waveform(t[Hz])');
%     model.physics('acpr').feature('pr1').set('p0', 1, 'p_sound*sin(2*pi*freq*t[Hz])');
%     model.physics('acpr').feature('pr1').set('p0', 1, 'p_sound*sin(2*pi*freq*t)');
    
end
%USEFUL COMMANDS TO REMEMBER

%%% show commands
% help mli

% get model parameter
hist_geometry.param.get('freq');

% show tags for notdes and subnodes in the model
mphmodel(new_geometry);

%%% gui for looking at model
% mphnavigator(new_geometry)