# Simplify making xlsx.workbooks

xlsx.tablefy <- function(
		workbook.table,
		sheet.name,
		workbook.name,
		distrib.table = FALSE,
		sort.it = TRUE,
		style = TRUE ,
		style.cols = "avg_log2FC",
		veryhigh.rule = 2,
		high.rule = 1.5,
		low.rule = -1.5,
		verylow.rule = -2,
		type = "expression",
		start.Col = 1,
		start.Row = 1,
		col.Names = TRUE
){
	if(distrib.table == FALSE){
		if(sort.it == TRUE){
			workbook.table <- workbook.table %>% 
				mutate(
					sort_col = abs(avg_log2FC)
				) %>%
				mutate(
					sort_group = case_when(
						avg_log2FC >= veryhigh.rule ~ 1,
						avg_log2FC <= verylow.rule ~ 2,
						avg_log2FC >= high.rule & avg_log2FC < veryhigh.rule ~ 3,
						avg_log2FC <= low.rule & avg_log2FC > verylow.rule ~ 4,
						TRUE ~ 5
					)
				) %>% 
				arrange(as.numeric(cluster), sort_group, -sort_col) %>% 
				dplyr::select(avg_log2FC, p_val_adj, cluster, gene) %>% 
				group_by(cluster) %>%
				group_split() 
		}else{
			workbook.table <- workbook.table %>% 
				group_by(cluster) %>% 
				dplyr::select(avg_log2FC, p_val_adj, cluster, gene) %>% 
				group_split()
		}
		
		row.max <- max(
			sapply(1:length(workbook.table), 
						 function(x){
						 	dim(workbook.table[[x]])[1]
						 }
			)
		)
		
		for(i in 1:length(workbook.table)){
			gene.count <- dim(workbook.table[[i]])[1]
			spacer <- row.max - gene.count
			if(spacer != 0){
				spacer <- rep(NA, spacer)
				workbook.table[[i]] <- workbook.table[[i]] %>%
					tibble::add_row(avg_log2FC = spacer, p_val_adj = spacer, cluster = spacer, gene = spacer, .after = gene.count)
			}
		}
		workbook.table <- workbook.table %>% 
			purrr::map(
				~ {
					cluster_name <- unique(.x$cluster)[1]
					rename_with(.x, ~ paste0(., "_", cluster_name))
				}
			) %>% 
			bind_cols() %>% 
			mutate_at(vars(grep("avg_log2FC|p_val_adj", colnames(workbook.table))), as.numeric)
		
		
		openxlsx::addWorksheet(workbook.name, sheetName = sheet.name)
		openxlsx::writeData(workbook.name, sheet = sheet.name, x = workbook.table, startCol = start.Col, startRow = start.Row, colNames = col.Names)
		
		if(style == TRUE){
			style.veryhigh <- openxlsx::createStyle(fontColour = "#FFFFFF", bgFill = "#377D43", textDecoration = "bold")
			style.high <- openxlsx::createStyle(fontColour = "#377D43", bgFill = "#CEEED0")
			style.low <- openxlsx::createStyle(fontColour = "#8E1C12", bgFill = "#F6C9CE")
			style.verylow <- openxlsx::createStyle(fontColour = "#FFFFFF", bgFill = "#8E1C12", textDecoration = "bold")
			
			if(is.character(style.cols)){
				style.cols <- grep(style.cols, colnames(workbook.table))
			}
			for(col.to.format in style.cols){
				openxlsx::conditionalFormatting(workbook.name, sheet = sheet.name, 
																				cols = col.to.format, 
																				rows = 1:nrow(workbook.table) +1,
																				rule =  c(high.rule, veryhigh.rule), 
																				style = style.high, 
																				type = "between", stack = TRUE)
				openxlsx::conditionalFormatting(workbook.name, sheet = sheet.name, 
																				cols = col.to.format, 
																				rows = 1:nrow(workbook.table) +1,
																				rule = paste0(">=", veryhigh.rule), 
																				style = style.veryhigh, 
																				type = "expression", stack = TRUE)
				openxlsx::conditionalFormatting(workbook.name, sheet = sheet.name, 
																				cols = col.to.format, 
																				rows = 1:nrow(workbook.table) +1,
																				rule = c(low.rule, verylow.rule), 
																				style = style.low, 
																				type = "between", stack = TRUE)
				openxlsx::conditionalFormatting(workbook.name, sheet = sheet.name, 
																				cols = col.to.format, 
																				rows = 1:nrow(workbook.table) +1,
																				rule = paste0("<=", verylow.rule), 
																				style = style.verylow, 
																				type = "expression", stack = TRUE)
			}
			
		} 
		
	} else if (distrib.table == TRUE){
		openxlsx::addWorksheet(
			workbook.name, 
			sheetName = sheet.name
		)
		openxlsx::writeData(
			workbook.name, 
			sheet = sheet.name, 
			x = workbook.table, 
			startCol = start.Col, startRow = start.Row, 
			colNames = col.Names
		)
		
	}
}
