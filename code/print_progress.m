function reverseStr = print_progress(current, total, reverseStr)
	percentDone = 100 * current / total;
    msg = sprintf('%3.2f percent done', percentDone);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
	if current >= total
		fprintf('\n')
	end


