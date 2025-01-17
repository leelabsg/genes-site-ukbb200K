'use strict';

$.getJSON('/static/genes.json').then(function(data) {
    $(document).ready(function() {
        var table = new Tabulator('#table', {
//            height: 600, // setting height lets Tabulator's VirtualDOM load really fast but makes scrolling awkward
            layout: 'fitColumns',
            pagination: 'local',
            paginationSize: 100,
            columns: [
                {title: 'Name', field:'name', formatter:'link', formatterParams: { urlPrefix: '/gene/' }, headerFilter:true, widthGrow:2},
                {title: 'Chrom', field:'chrom', formatter:'link', headerFilter:true, headerFilterFunc:'='},
                {title: 'Num p&le;2.5<sup>-6</sup> Associations', field:'num_sig_assocs', formatter:'comma_fmt', width:200},
                {title: 'Best Phenotype', field:'best_assoc'},
                {title: 'Best P-value', field:'best_pval', formatter:'2digit_fmt'},
            ],
            data: data,
            initialSort: [{column:"num_sig_assocs", dir:"desc"}],
            tooltipGenerationMode:'hover',tooltips:tabulator_tooltip_maker,tooltipsHeader:true,
        });
    });
});
