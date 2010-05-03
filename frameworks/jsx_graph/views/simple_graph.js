// ==========================================================================
// Project:   JsxGraph.SimpleGraphView
// Copyright: ©2010 My Company, Inc.
// ==========================================================================
/*globals JsxGraph */

/** @class

  (Document Your View Here)

  @extends SC.View
*/
JsxGraph.SimpleGraphView = SC.View.extend(
/** @scope JsxGraph.SimpleGraphView.prototype */ {

  // TODO: Add your own code here.
  
  classNames: "jxgbox",
  
  layerDidChange: function() {
    this.set('layerNeedsUpdate', YES);
  }.observes('layer'),

  updateLayer: function() {
    sc_super();
    var layer = this.get('layer');
    if (layer) {
      JXG.JSXGraph.initBoard(layer.id, {originX: 250, originY: 250, unitX: 50, unitY: 50, axis:true});
    }
  }

});
