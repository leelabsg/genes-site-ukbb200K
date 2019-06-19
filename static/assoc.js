'use strict';

$.getJSON('/api/variants/'+model.genename+'/'+model.phecode).then(function(resp) {
    window._debug.resp = resp;
    var df = resp.df;
    df.pos = df.pos_delta;
    delete df.pos_delta;
    _.map(df.pos, function(elem, i) { if (i>0) { df.pos[i] += df.pos[i-1]; } })
    var data = dataframe_to_objects(df);
    console.log(data);

    $(function() {
        var table = new Tabulator('#table', {
            //height: 600, // setting height lets Tabulator's VirtualDOM load really fast but makes scrolling awkward
            layout: 'fitColumns',
            pagination: 'local', // TODO: disable pagination if <100 variants
            paginationSize: 100,
            columns: [
                {title: 'Position on chr'+resp.chrom, field:'pos'},
                {title: 'Allele', field:'base'},
                {title: 'MAF (Minor Allele Frequency)', field:'maf'},
                {title: 'Case MAC (Minor Allele Count)', field:'mac_case'},
                {title: 'Control MAC (Minor Allele Count)', field:'mac_control'},
                {title: 'P-value', field:'pval'},
            ],
            data: data,
            initialSort: [{column:'pval', dir:'asc'}],
            tooltipGenerationMode:'hover',tooltips:tabulator_tooltip_maker,tooltipsHeader:true,
        });
    });
});
