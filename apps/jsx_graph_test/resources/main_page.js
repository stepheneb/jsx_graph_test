// ==========================================================================
// Project:   JsxGraphTest - mainPage
// Copyright: ©2010 My Company, Inc.
// ==========================================================================
/*globals JsxGraphTest */

// This page describes the main user interface for your application.  
JsxGraphTest.mainPage = SC.Page.design({

  // The main pane is made visible on screen as soon as your app is loaded.
  // Add childViews to this pane for views to display immediately on page 
  // load.
  mainPane: SC.MainPane.design({
    childViews: 'tabView'.w(),
    
    tabView: SC.TabView.design({ 
      layout: {top: 30, bottom: 5, left: 5, right: 5 }, 
      items: [ 
        {title: "Welcome", value: "demos" },
        {title: "JSXGraph", value: "jsxGraphViews" },
      ], 
      itemTitleKey: 'title', 
      itemValueKey: 'value', 
      nowShowing: 'demos', // defining the startup tab 
      userDefaultKey: 'mainPaneTab'
    })

  }),
  
  demos: SC.LabelView.design({
    escapeHTML: NO,
    value: "<h1>Welcome</h1><p>Some initial experiments embedding <a href='http://jsxgraph.uni-bayreuth.de/wp/examples/'>JSXGraph</a> into <a href='http://www.sproutcore.com/'>SproutCore</a>.</p>",
    layout: { centerX: 0, centerY: 0, width: 500, height: 400 }
  }),
  
  jsxGraphViews: SC.View.design({
    childViews: 'jsxGraphTabsView'.w(),
    
    jsxGraphTabsView: SC.TabView.design({ 
      layout: {top: 30, bottom: 5, left: 5, right: 5 }, 
      items: [ 
        {title: "Empty Graph", value: "jsxGraph1" },
        {title: "Curves", value: "jsxGraph2" },
        {title: "Bezier Curves", value: "jsxGraph3" },
        {title: "SimpleGraph", value: "jsxGraph4" }
      ], 
      itemTitleKey: 'title', 
      itemValueKey: 'value', 
      nowShowing: 'jsxGraph1', // defining the startup tab 
      userDefaultKey: 'jsxGraphTab'
    })
  }),

  jsxGraph1: SC.View.design({
    layout: {width:500, height:500, centerX:0, centerY:0},
    classNames: "jxgbox",
    layerDidChange: function() {
      this.set('layerNeedsUpdate', YES);
    }.observes('layer'),
  
    updateLayer: function() {
      sc_super();
      layer = this.get('layer');
      if (layer) {
        JXG.JSXGraph.initBoard(layer.id, {originX: 250, originY: 250, unitX: 50, unitY: 50, axis:true});
      }
    }
  }),

  jsxGraph2: SC.View.design({
    layout: {width:500, height:500, centerX:0, centerY:0},
    classNames: "jxgbox",
    layerDidChange: function() {
      this.set('layerNeedsUpdate', YES);
    }.observes('layer'),
  
    updateLayer: function() {
      sc_super();
      layer = this.get('layer');
      if (layer) {
        var brd = JXG.JSXGraph.initBoard(layer.id,
                   {axis:true,originX: 250, originY: 250, unitX: 50, unitY: 25});
        brd.suspendUpdate();
        var p = [];
        p[0] = brd.create('point', [-4,2], {style:6});
        p[1] = brd.create('point', [3,-1], {style:6});
        p.push(brd.create('point', [-2,(Math.random()-0.5)*3],{style:6}));
        p.push(brd.create('point', [0.5,(Math.random()-0.5)*3],{style:6}));
        p.push(brd.create('point', [1,(Math.random()-0.5)*3],{style:6}));
        brd.update();
        var pol = brd.lagrangePolynomial(p);
        var g = brd.create('functiongraph', [pol, -10, 10], {strokeWidth:3});
        var g2 = brd.create('functiongraph', [brd.D(pol), -10, 10],
                                              {dash:3, strokeColor:'#ff0000'});
        brd.unsuspendUpdate();
      }
    }
  }),
  
  jsxGraph3: SC.View.design({
    layout: {width:500, height:500, centerX:0, centerY:0},
    classNames: "jxgbox",
    layerDidChange: function() {
      this.set('layerNeedsUpdate', YES);
    }.observes('layer'),
  
    updateLayer: function() {
      sc_super();
      layer = this.get('layer');
      if (layer) {
        
        var brd = JXG.JSXGraph.initBoard(layer.id,{boundingbox:[-4,4,4,-4],keepaspectratio:true,axis:true});

        var p = [];

        var col = 'red'; 
        p.push(brd.create('point',[2,1],{strokeColor:col,fillColor:col}));        // data point
        col = 'blue'; 
        p.push(brd.create('point',[0.75,2.5],{strokeColor:col,fillColor:col}));   // control point
        p.push(brd.create('point',[-0.3,0.3],{strokeColor:col,fillColor:col}));   // control point

        col = 'red'; 
        p.push(brd.create('point',[-3,1],{strokeColor:col,fillColor:col}));       // data point
        col = 'blue'; 
        p.push(brd.create('point',[-0.75,-2.5],{strokeColor:col,fillColor:col})); // control point
        p.push(brd.create('point',[1.5,-2.8],{strokeColor:col,fillColor:col}));   // control point

        col = 'red'; 
        p.push(brd.create('point',[2,-0.5],{strokeColor:col,fillColor:col}));     // data point

        var c = brd.create('curve', JXG.Math.Numerics.bezier(p), 
                       {strokecolor:'blue', strokeOpacity:0.6, strokeWidth:5}); 
      }
    }
  }),

  jsxGraph4: JsxGraph.SimpleGraphView.design({
    layout: {width:500, height:500, centerX:0, centerY:0}
  })

});
