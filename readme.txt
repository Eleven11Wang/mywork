class brainspanwork 
    class sample_charactor(brainspanwork)
        function : filter the tissue that less then five 
        function : find which columns is which tissue (find_columns_tissue)
        function : seperate the matrix by age (seperate_by_time)
        seperate : the matrix by sex (seperate_by_sex)
    (not done)
    class expression_deal_with(sample_charactor)
        function : tissue expression matrix  return a work matrix(tissue_matrix_workon)

        function : 50% of the sample have this gene return a filted dataframe 
        function : check normalized or not ,box plot of data
        function : upper quantile transformation 
        function : log transform 

    class expression_analysis(expression_deal_with)        
        function : differental expression age and sex 
        function : pca 
        function : clustering 
        function : pearson correlation 

    class network_analysis(expression_analysis):
        function : on class / paper research 

