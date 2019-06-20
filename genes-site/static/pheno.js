'use strict';

$.getJSON('/api/pheno/'+model.phecode).then(function(resp) {
    $(function() {
        var unbinned_variants = dataframe_to_objects(resp.assocs);
        var variant_bins = resp.manhattan_bins;
        unbinned_variants.forEach(function(v){v.pos=v.startpos});
        var significance_threshold = 0.05 / resp.num_genes;
        console.log(unbinned_variants);
        console.log(variant_bins);
        create_gwas_plot(variant_bins, unbinned_variants, significance_threshold);
    });

    var num_genes = resp.num_genes;
    var assocs = objects_to_dataframe(_.sortBy(dataframe_to_objects(resp.assocs), _.property('startpos')));

    assocs.id = assocs.name;
    assocs.trait_label = assocs.name;
    assocs.log_pvalue = assocs.pval.map(function(p) { return -Math.log10(Math.max(1e-100, p)); });
    assocs.trait_group = assocs.chrom.map(function(s) { return 'chr'+s.toString().padStart(2,'0'); });

    var significance_threshold = -Math.log10(0.05 / num_genes);
    var best_nlpval = d3.max(assocs.log_pvalue);

    var data_sources = new LocusZoom.DataSources().add('phewas', ['StaticJSON', assocs]);
    var layout = {
        width: 800,
        height: 400,
        min_width: 800,
        min_height: 400,
        responsive_resize: 'width_only',
        mouse_guide: false,
        dashboard: {components: [ {type: 'download', position: 'right', color: 'gray' } ]},
        panels: [
            LocusZoom.Layouts.get('panel', 'phewas', {
                margin: {top: 5, right: 5, bottom: 50, left: 50 }
            })
        ],
    }

    layout.panels[0].data_layers[0].offset = significance_threshold;
    layout.panels[0].data_layers[1].fields.push('phewas:num_rare');
    layout.panels[0].data_layers[1].fields.push('phewas:chrom', 'phewas:startpos', 'phewas:endpos');
    layout.panels[0].data_layers[1].fields.push('phewas:mac_case', 'phewas:mac_control');
    layout.panels[0].data_layers[1].tooltip.html =
        ("gene: <strong>{{phewas:trait_label|htmlescape}}</strong><br>" +
         "P-value: <strong>{{phewas:log_pvalue|logtoscinotation|htmlescape}}</strong><br>" +
         "#Rare Variants: <strong>{{phewas:num_rare}}</strong><br>" +
         "Case / Control MAC: <strong>{{phewas:mac_case}} / {{phewas:mac_control}}</strong><br>" +
         "Chrom:Start-End: <strong>{{phewas:chrom}} : {{phewas:startpos}} - {{phewas:endpos}}</strong><br>"
        );
    layout.panels[0].data_layers[1].behaviors.onclick = [{action: 'link', href: '/assoc/{{phewas:id}}/'+model.phecode}];
    layout.panels[0].data_layers[1].y_axis.min_extent = [0, significance_threshold*1.1];

    if (assocs.id.length <= 10) {
        layout.panels[0].data_layers[1].label.filters = []; // show all labels
    } else if (assocs.log_pvalue.filter(function(nlpval) { return nlpval == best_nlpval; }).length >= 6) {
        layout.panels[0].data_layers[1].label = false; // too many are tied for 1st and will make a mess so just hide all labels
    } else {
        var eighth_best_nlpval = _.sortBy(assocs.log_pvalue).reverse()[8];
        layout.panels[0].data_layers[1].label.filters = [
            {field: 'phewas:log_pvalue', operator: '>', value: significance_threshold},
            {field: 'phewas:log_pvalue', operator: '>', value: best_nlpval*0.5}, // must be in top half of screen
            {field: 'phewas:log_pvalue', operator: '>=', value: eighth_best_nlpval} // must be among the best
        ];
    }

    window._debug.assocs = assocs;
    $(function() {
        var plot = LocusZoom.populate("#phewas_plot_container", data_sources, layout);
        window._debug.plot = plot;
    });

    $(function() {
        var data = dataframe_to_objects(assocs);
        var table = new Tabulator('#table', {
            //height: 600, // setting height lets Tabulator's VirtualDOM load really fast but makes scrolling awkward
            layout: 'fitColumns',
            pagination: 'local',
            paginationSize: 15,
            columns: [
                {title: 'Gene', field:'name', formatter:'link', formatterParams: {url:function(cell){return '/assoc/'+cell.getValue()+'/'+model.phecode}}, headerFilter:true, widthGrow:3},
                {title: 'P-value', field:'pval'},
                {title: '#Rare Variants', field:'num_rare'},
                {title: 'Chromosome', field:'chrom', headerFilter:true, headerFilterFunc:'='},
                {title: 'Start-End', field:'startpos', formatter: function(cell){return cell.getValue()+' - '+cell.getData().endpos}},
                {title: 'Case MAC (Minor Allele Count)', field:'mac_case'},
                {title: 'Control MAC (Minor Allele Count)', field:'mac_control'},
            ],
            data: data,
            initialSort: [{column:'pval', dir:'asc'}],
            tooltipGenerationMode:'hover',tooltips:tabulator_tooltip_maker,tooltipsHeader:true,
        });
    });
});


function create_gwas_plot(variant_bins, unbinned_variants, significance_threshold) { // from PheWeb

    if (variant_bins.length && typeof variant_bins[0].qvals === "undefined") {
        // this json was generated by an old version of pheweb, so we'll manually fix things up.
        variant_bins.forEach(function(bin) {
            bin.qvals = bin.neglog10_pvals || bin.nlpvals; delete bin.neglog10_pvals; delete bin.nlpvals;
            bin.qval_extents = bin.neglog10_pval_extents || bin.nlpval_extents; delete bin.neglog10_pval_extents; delete bin.nlpval_extents;
        });
    }

    var get_chrom_offsets = _.memoize(function() {
        var chrom_padding = 2e7;
        var chrom_extents = {};

        var update_chrom_extents = function(variant) {
            if (!(variant.chrom in chrom_extents)) {
                chrom_extents[variant.chrom] = [variant.pos, variant.pos];
            } else if (variant.pos > chrom_extents[variant.chrom][1]) {
                chrom_extents[variant.chrom][1] = variant.pos;
            } else if (variant.pos < chrom_extents[variant.chrom][0]) {
                chrom_extents[variant.chrom][0] = variant.pos;
            }
        }
        variant_bins.forEach(update_chrom_extents);
        unbinned_variants.forEach(update_chrom_extents);

        var chroms = _.sortBy(Object.keys(chrom_extents), parseInt);

        var chrom_genomic_start_positions = {};
        chrom_genomic_start_positions[chroms[0]] = 0;
        for (var i=1; i<chroms.length; i++) {
            chrom_genomic_start_positions[chroms[i]] = chrom_genomic_start_positions[chroms[i-1]] + chrom_extents[chroms[i-1]][1] - chrom_extents[chroms[i-1]][0] + chrom_padding;
        }

        // chrom_offsets are defined to be the numbers that make `get_genomic_position()` work.
        // ie, they leave a gap of `chrom_padding` between the last variant on one chromosome and the first on the next.
        var chrom_offsets = {};
        Object.keys(chrom_genomic_start_positions).forEach(function(chrom) {
            chrom_offsets[chrom] = chrom_genomic_start_positions[chrom] - chrom_extents[chrom][0];
        });

        return {
            chrom_extents: chrom_extents,
            chroms: chroms,
            chrom_genomic_start_positions: chrom_genomic_start_positions,
            chrom_offsets: chrom_offsets,
        };
    });

    function get_genomic_position(variant) {
        var chrom_offsets = get_chrom_offsets().chrom_offsets;
        return chrom_offsets[variant.chrom] + variant.pos;
    }

    function get_y_axis_config(max_data_qval, plot_height, includes_pval0) {

        var possible_ticks = [];
        if (max_data_qval <= 14) { possible_ticks = _.range(0, 14.1, 2); }
        else if (max_data_qval <= 28) { possible_ticks = _.range(0, 28.1, 4); }
        else if (max_data_qval <= 40) { possible_ticks = _.range(0, 40.1, 8); }
        else {
            possible_ticks = _.range(0, 20.1, 4);
            if (max_data_qval <= 70) { possible_ticks = possible_ticks.concat([30,40,50,60,70]); }
            else if (max_data_qval <= 120) { possible_ticks = possible_ticks.concat([40,60,80,100,120]); }
            else if (max_data_qval <= 220) { possible_ticks = possible_ticks.concat([60,100,140,180,220]); }
            else {
                var power_of_ten = Math.pow(10, Math.floor(Math.log10(max_data_qval)));
                var first_digit = max_data_qval / power_of_ten;
                var multipliers;
                if (first_digit <= 2) { multipliers = [0.5, 1, 1.5, 2]; }
                else if (first_digit <= 4) { multipliers = [1, 2, 3, 4]; }
                else { multipliers = [2, 4, 6, 8, 10]; }
                possible_ticks = possible_ticks.concat(multipliers.map(function(m) { return m * power_of_ten; }));
            }
        }
        // Include all ticks < qval.  Then also include the next tick.
        // That should mean we'll always have the largest tick >= the largest variant.
        var ticks = possible_ticks.filter(function(qval) { return qval < max_data_qval; });
        if (ticks.length < possible_ticks.length) { ticks.push(possible_ticks[ticks.length]); }

        // Use the largest tick for the top of our y-axis so that we'll have a tick nicely rendered right at the top.
        var max_plot_qval = ticks[ticks.length-1];
        // If we have any qval=inf (pval=0) variants, leave space for them.
        if (includes_pval0) { max_plot_qval *= 1.1 }
        var scale = d3.scale.linear().clamp(true);
        if (max_plot_qval <= 40) {
            scale = scale
                .domain([max_plot_qval, 0])
                .range([0, plot_height]);
        } else {
            scale = scale
                .domain([max_plot_qval, 20, 0])
                .range([0, plot_height/2, plot_height]);
        }

        if (includes_pval0) { ticks.push(Infinity); }

        return {
            'scale': scale,
            'draw_break_at_20': !(max_plot_qval <= 40),
            'ticks': ticks,
        };
    }

    $(function() {
        var svg_width = $('#manhattan_plot_container').width();
        var svg_height = 550;
        var plot_margin = {
            'left': 70,
            'right': 30,
            'top': 20,
            'bottom': 50,
        };
        var plot_width = svg_width - plot_margin.left - plot_margin.right;
        var plot_height = svg_height - plot_margin.top - plot_margin.bottom;

        var gwas_svg = d3.select('#manhattan_plot_container').append("svg")
            .attr('id', 'gwas_svg')
            .attr("width", svg_width)
            .attr("height", svg_height)
            .style("display", "block")
            .style("margin", "auto");
        var gwas_plot = gwas_svg.append("g")
            .attr('id', 'gwas_plot')
            .attr("transform", fmt("translate({0},{1})", plot_margin.left, plot_margin.top));

        // Significance Threshold line
        var significance_threshold_tooltip = d3.tip()
            .attr('class', 'd3-tip')
            .html('Significance Threshold: '+significance_threshold.toExponential(1))
            .offset([-8,0]);
        gwas_svg.call(significance_threshold_tooltip);

        var genomic_position_extent = (function() {
            var extent1 = d3.extent(variant_bins, get_genomic_position);
            var extent2 = d3.extent(unbinned_variants, get_genomic_position);
            return d3.extent(extent1.concat(extent2));
        })();

        var x_scale = d3.scale.linear()
            .domain(genomic_position_extent)
            .range([0, plot_width]);

        var includes_pval0 = _.any(unbinned_variants, function(variant) { return variant.pval === 0; });

        var highest_plot_qval = Math.max(
            -Math.log10(significance_threshold) + 0.5,
            (function() {
                var best_unbinned_qval = -Math.log10(d3.min(unbinned_variants, function(d) {
                    return (d.pval === 0) ? 1 : d.pval;
                }));
                if (best_unbinned_qval !== undefined) return best_unbinned_qval;
                return d3.max(variant_bins, function(bin) {
                    return d3.max(bin, _.property('qval'));
                });
            })());

        var y_axis_config = get_y_axis_config(highest_plot_qval, plot_height, includes_pval0);
        var y_scale = y_axis_config.scale;

        // TODO: draw a small y-axis-break at 20 if `y_axis_config.draw_break_at_20`
        var y_axis = d3.svg.axis()
            .scale(y_scale)
            .orient("left")
            .tickFormat(d3.format("d"))
            .tickValues(y_axis_config.ticks)
        gwas_plot.append("g")
            .attr("class", "y axis")
            .attr('transform', 'translate(-8,0)') // avoid letting points spill through the y axis.
            .call(y_axis);

        if (includes_pval0) {
            var y_axis_break_inf_offset = y_scale(Infinity) + (y_scale(0)-y_scale(Infinity)) * 0.03
            gwas_plot.append('line')
                .attr('x1', -8-7).attr('x2', -8+7)
                .attr('y1', y_axis_break_inf_offset+6).attr('y2', y_axis_break_inf_offset-6)
                .attr('stroke', '#666').attr('stroke-width', '3px');
        }
        if (y_axis_config.draw_break_at_20) {
            var y_axis_break_20_offset = y_scale(20);
            gwas_plot.append('line')
                .attr('x1', -8-7).attr('x2', -8+7)
                .attr('y1', y_axis_break_20_offset+6).attr('y2', y_axis_break_20_offset-6)
                .attr('stroke', '#666').attr('stroke-width', '3px');
        }

        gwas_svg.append('text')
            .style('text-anchor', 'middle')
            .attr('transform', fmt('translate({0},{1})rotate(-90)',
                                   plot_margin.left*.4,
                                   plot_height/2 + plot_margin.top))
            .text('-log\u2081\u2080(p-value)'); // Unicode subscript "10"

        var chroms_and_midpoints = (function() {
            var v = get_chrom_offsets();
            return v.chroms.map(function(chrom) {
                return {
                    chrom: chrom,
                    midpoint: v.chrom_genomic_start_positions[chrom] + (v.chrom_extents[chrom][1] - v.chrom_extents[chrom][0]) / 2,
                };
            });
        })();

        var color_by_chrom = d3.scale.ordinal()
            .domain(get_chrom_offsets().chroms)
            .range(['rgb(120,120,186)', 'rgb(0,0,66)']);
        //colors to maybe sample from later:
        //.range(['rgb(120,120,186)', 'rgb(0,0,66)', 'rgb(44,150,220)', 'rgb(40,60,80)', 'rgb(33,127,188)', 'rgb(143,76,176)']);

        gwas_svg.selectAll('text.chrom_label')
            .data(chroms_and_midpoints)
            .enter()
            .append('text')
            .style('text-anchor', 'middle')
            .attr('transform', function(d) {
                return fmt('translate({0},{1})',
                           plot_margin.left + x_scale(d.midpoint),
                           plot_height + plot_margin.top + 20);
            })
            .text(function(d) {
                return d.chrom;
            })
            .style('fill', function(d) {
                return color_by_chrom(d.chrom);
            });

        gwas_plot.append('line')
            .attr('x1', 0)
            .attr('x2', plot_width)
            .attr('y1', y_scale(-Math.log10(significance_threshold)))
            .attr('y2', y_scale(-Math.log10(significance_threshold)))
            .attr('stroke-width', '5px')
            .attr('stroke', 'lightgray')
            .attr('stroke-dasharray', '10,10')
            .on('mouseover', significance_threshold_tooltip.show)
            .on('mouseout', significance_threshold_tooltip.hide);

        // Points & labels
        var tooltip_template = _.template(
            '<%= d.name %><br>' +
                'P-value: <%= d.pval %><br>' +
                '#Rare Variants: <%= d.num_rare %><br>' +
                'Case / Control MAC: <%= d.mac_case %> / <%= d.mac_control %><br>' +
                'Chrom:Start-End: <%= d.chrom %>:<%= d.startpos %>-<%= d.endpos %>' +
                '');
        var point_tooltip = d3.tip()
            .attr('class', 'd3-tip')
            .html(function(d) {
                return tooltip_template({d: d});
            })
            .offset([-6,0]);
        gwas_svg.call(point_tooltip);

        function get_link_to_LZ(d) {
            return fmt('/assoc/{0}/{1}', d.name, model.phecode);
        }

        // TODO: if the label touches any circles or labels, skip it?
        var variants_to_label = _.sortBy(unbinned_variants, _.property('pval'))
            .filter(function(d) { return d.pval < significance_threshold; })
            .slice(0,7);
        var genenames = gwas_plot.append('g')
            .attr('class', 'genenames')
            .selectAll('text.genenames')
            .data(variants_to_label)
            .enter()
            .append('text')
            .attr('class', 'genename_text')
            .style('font-style', 'italic')
            .attr('text-anchor', 'middle')
            .attr('transform', function(d) {
                return fmt('translate({0},{1})',
                           x_scale(get_genomic_position(d)),
                           y_scale(-Math.log10(d.pval))-5);
            })
            .text(function(d) {
                return d.name;
            });

        function pp1() {
        gwas_plot.append('g')
            .attr('class', 'variant_hover_rings')
            .selectAll('a.variant_hover_ring')
            .data(unbinned_variants)
            .enter()
            .append('a')
            .attr('class', 'variant_hover_ring')
            .attr('xlink:href', get_link_to_LZ)
            .append('circle')
            .attr('cx', function(d) {
                return x_scale(get_genomic_position(d));
            })
            .attr('cy', function(d) {
                return y_scale(-Math.log10(d.pval));
            })
            .attr('r', 7)
            .style('opacity', 0)
            .style('stroke-width', 1)
            .on('mouseover', function(d) {
                //Note: once a tooltip has been explicitly placed once, it must be explicitly placed forever after.
                var target_node = document.getElementById(fmt('variant-point-{0}-{1}-{2}-{3}', d.chrom, d.pos, d.ref, d.alt));
                point_tooltip.show(d, target_node);
            })
            .on('mouseout', point_tooltip.hide);
        }
        pp1();

        function pp2() {
        gwas_plot.append('g')
            .attr('class', 'variant_points')
            .selectAll('a.variant_point')
            .data(unbinned_variants)
            .enter()
            .append('a')
            .attr('class', 'variant_point')
            .attr('xlink:href', get_link_to_LZ)
            .append('circle')
            .attr('id', function(d) {
                return fmt('variant-point-{0}-{1}-{2}-{3}', d.chrom, d.pos, d.ref, d.alt);
            })
            .attr('cx', function(d) {
                return x_scale(get_genomic_position(d));
            })
            .attr('cy', function(d) {
                return y_scale(-Math.log10(d.pval));
            })
            .attr('r', 2.3)
            .style('fill', function(d) {
                return color_by_chrom(d.chrom);
            })
            .on('mouseover', function(d) {
                //Note: once a tooltip has been explicitly placed once, it must be explicitly placed forever after.
                point_tooltip.show(d, this);
            })
            .on('mouseout', point_tooltip.hide);
        }
        pp2();

        function pp3() { // drawing the ~60k binned variant circles takes ~500ms.  The (far fewer) unbinned variants take much less time.
        var bins = gwas_plot.append('g')
            .attr('class', 'bins')
            .selectAll('g.bin')
            .data(variant_bins)
            .enter()
            .append('g')
            .attr('class', 'bin')
            .each(function(d) { //todo: do this in a forEach
                d.x = x_scale(get_genomic_position(d));
                d.color = color_by_chrom(d.chrom);
            });
        bins.selectAll('circle.binned_variant_point')
            .data(_.property('qvals'))
            .enter()
            .append('circle')
            .attr('class', 'binned_variant_point')
            .attr('cx', function(d, i, parent_i) {
                return variant_bins[parent_i].x;
            })
            .attr('cy', function(qval) {
                return y_scale(qval);
            })
            .attr('r', 2.3)
            .style('fill', function(d, i, parent_i) {
                // return color_by_chrom(d3.select(this.parentNode).datum().chrom); //slow
                // return color_by_chrom(this.parentNode.__data__.chrom); //slow?
                // return this.parentNode.__data__.color;
                return variant_bins[parent_i].color;
            });
        bins.selectAll('circle.binned_variant_line')
            .data(_.property('qval_extents'))
            .enter()
            .append('line')
            .attr('class', 'binned_variant_line')
            .attr('x1', function(d, i, parent_i) { return variant_bins[parent_i].x; })
            .attr('x2', function(d, i, parent_i) { return variant_bins[parent_i].x; })
            .attr('y1', function(d) { return y_scale(d[0]); })
            .attr('y2', function(d) { return y_scale(d[1]); })
            .style('stroke', function(d, i, parent_i) { return variant_bins[parent_i].color; })
            .style('stroke-width', 4.6)
            .style('stroke-linecap', 'round');
        }
        pp3();

    });
}
