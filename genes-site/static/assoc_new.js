'use strict';

//window.alert(model.associd)
$.getJSON('/api/assocgroup/'+model.associd).then(function(resp) {
//$.getJSON('/api/assocgroup/1').then(function(resp) {
    
	var data = dataframe_to_objects(resp.assocgroup);
	$(document).ready(function() {
        var table = new Tabulator('#table', {
//            height: 600, // setting height lets Tabulator's VirtualDOM load really fast but makes scrolling awkward
            layout: 'fitColumns',
            pagination: 'local',
            paginationSize: 9,
            columns: [
                {title: 'Group', field:'description'},
                {title: 'P-value', field:'pval', formatter:'2digit_fmt'},
            ],
			data: data,
			initialSort: [{column:'description', dir:'asc'}],
            tooltipGenerationMode:'hover',tooltips:tabulator_tooltip_maker,tooltipsHeader:true,
        });
    });
});